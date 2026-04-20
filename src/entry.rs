use crate::constants::{MODE_AREA, MODE_P2P};
use crate::diffraction::diffraction_loss;
use crate::enums::*;
use crate::errors::*;
use crate::helper::{fortran_dim, itm_max, itm_min};
use crate::initialization::{initialize_area, initialize_point_to_point, quick_pfl};
use crate::los::line_of_sight_loss;
use crate::troposcatter::troposcatter_loss;
use crate::types::{ComplexDouble, IntermediateValues};
use crate::variability::{free_space_loss, validate_inputs, variability};

pub fn longley_rice(
    theta_hzn: [f64; 2],
    f__mhz: f64,
    z_g: ComplexDouble,
    d_hzn__meter: [f64; 2],
    h_e__meter: [f64; 2],
    gamma_e: f64,
    n_s: f64,
    delta_h__meter: f64,
    h__meter: [f64; 2],
    d__meter: f64,
    mode: i32,
    a_ref__db: &mut f64,
    warnings: &mut i32,
    propmode: &mut i32,
) -> i32 {
    let a_e__meter = 1.0 / gamma_e;

    let mut d_hzn_s__meter = [0.0f64; 2];
    for i in 0..2 {
        d_hzn_s__meter[i] = (2.0 * h_e__meter[i] * a_e__meter).sqrt();
    }

    let d_sML__meter = d_hzn_s__meter[0] + d_hzn_s__meter[1];
    let d_ML__meter = d_hzn__meter[0] + d_hzn__meter[1];

    let theta_los = -itm_max(theta_hzn[0] + theta_hzn[1], -d_ML__meter / a_e__meter);

    if theta_hzn[0].abs() > 200e-3 {
        *warnings |= WARN_TX_HORIZON_ANGLE;
    }
    if theta_hzn[1].abs() > 200e-3 {
        *warnings |= WARN_RX_HORIZON_ANGLE;
    }

    if d_hzn__meter[0] < 0.1 * d_hzn_s__meter[0] {
        *warnings |= WARN_TX_HORIZON_DISTANCE_1;
    }
    if d_hzn__meter[1] < 0.1 * d_hzn_s__meter[1] {
        *warnings |= WARN_RX_HORIZON_DISTANCE_1;
    }

    if d_hzn__meter[0] > 3.0 * d_hzn_s__meter[0] {
        *warnings |= WARN_TX_HORIZON_DISTANCE_2;
    }
    if d_hzn__meter[1] > 3.0 * d_hzn_s__meter[1] {
        *warnings |= WARN_RX_HORIZON_DISTANCE_2;
    }

    if n_s < 150.0 {
        return ERROR_SURFACE_REFRACTIVITY_SMALL;
    }
    if n_s > 400.0 {
        return ERROR_SURFACE_REFRACTIVITY_LARGE;
    }
    if n_s < 250.0 {
        *warnings |= WARN_SURFACE_REFRACTIVITY;
    }

    if !(4000000.0..=13333333.0).contains(&a_e__meter) {
        return ERROR_EFFECTIVE_EARTH;
    }

    if z_g.re <= z_g.im.abs() {
        return ERROR_GROUND_IMPEDANCE;
    }

    let d_3__meter = itm_max(
        d_sML__meter,
        d_ML__meter + 5.0 * (a_e__meter.powi(2) / f__mhz).powf(1.0 / 3.0),
    );
    let d_4__meter = d_3__meter + 10.0 * (a_e__meter.powi(2) / f__mhz).powf(1.0 / 3.0);

    let a_3__db = diffraction_loss(
        d_3__meter,
        d_hzn__meter,
        h_e__meter,
        z_g,
        a_e__meter,
        delta_h__meter,
        h__meter,
        mode,
        theta_los,
        d_sML__meter,
        f__mhz,
    );
    let a_4__db = diffraction_loss(
        d_4__meter,
        d_hzn__meter,
        h_e__meter,
        z_g,
        a_e__meter,
        delta_h__meter,
        h__meter,
        mode,
        theta_los,
        d_sML__meter,
        f__mhz,
    );

    let m_d = (a_4__db - a_3__db) / (d_4__meter - d_3__meter);
    let a_d0__db = a_3__db - m_d * d_3__meter;

    let d_min__meter = (h_e__meter[0] - h_e__meter[1]).abs() / 200e-3;

    if d__meter < d_min__meter {
        *warnings |= WARN_PATH_DISTANCE_TOO_SMALL_1;
    }
    if d__meter < 1e3 {
        *warnings |= WARN_PATH_DISTANCE_TOO_SMALL_2;
    }
    if d__meter > 1000e3 {
        *warnings |= WARN_PATH_DISTANCE_TOO_BIG_1;
    }
    if d__meter > 2000e3 {
        *warnings |= WARN_PATH_DISTANCE_TOO_BIG_2;
    }

    if d__meter < d_sML__meter {
        let a_sML__db = d_sML__meter * m_d + a_d0__db;

        let mut d_0__meter = 0.04 * f__mhz * h_e__meter[0] * h_e__meter[1];

        let d_1__meter;
        if a_d0__db >= 0.0 {
            d_0__meter = itm_min(d_0__meter, 0.5 * d_ML__meter);
            d_1__meter = d_0__meter + 0.25 * (d_ML__meter - d_0__meter);
        } else {
            d_1__meter = itm_max(-a_d0__db / m_d, 0.25 * d_ML__meter);
        }

        let a_1__db = line_of_sight_loss(
            d_1__meter,
            h_e__meter,
            z_g,
            delta_h__meter,
            m_d,
            a_d0__db,
            d_sML__meter,
            f__mhz,
        );

        let mut flag = false;

        let mut k_hat_1__db_per_meter = 0.0f64;
        let mut k_hat_2__db_per_meter = 0.0f64;

        if d_0__meter < d_1__meter {
            let a_0__db = line_of_sight_loss(
                d_0__meter,
                h_e__meter,
                z_g,
                delta_h__meter,
                m_d,
                a_d0__db,
                d_sML__meter,
                f__mhz,
            );

            let q = (d_sML__meter / d_0__meter).ln();

            k_hat_2__db_per_meter = itm_max(
                0.0,
                ((d_sML__meter - d_0__meter) * (a_1__db - a_0__db)
                    - (d_1__meter - d_0__meter) * (a_sML__db - a_0__db))
                    / ((d_sML__meter - d_0__meter) * (d_1__meter / d_0__meter).ln()
                        - (d_1__meter - d_0__meter) * q),
            );

            flag = a_d0__db > 0.0 || k_hat_2__db_per_meter > 0.0;

            if flag {
                k_hat_1__db_per_meter =
                    (a_sML__db - a_0__db - k_hat_2__db_per_meter * q) / (d_sML__meter - d_0__meter);

                if k_hat_1__db_per_meter < 0.0 {
                    k_hat_1__db_per_meter = 0.0;
                    k_hat_2__db_per_meter = fortran_dim(a_sML__db, a_0__db) / q;

                    if k_hat_2__db_per_meter == 0.0 {
                        k_hat_1__db_per_meter = m_d;
                    }
                }
            }
        }

        if !flag {
            k_hat_1__db_per_meter = fortran_dim(a_sML__db, a_1__db) / (d_sML__meter - d_1__meter);
            k_hat_2__db_per_meter = 0.0;

            if k_hat_1__db_per_meter == 0.0 {
                k_hat_1__db_per_meter = m_d;
            }
        }

        let a_o__db = a_sML__db
            - k_hat_1__db_per_meter * d_sML__meter
            - k_hat_2__db_per_meter * d_sML__meter.ln();

        *a_ref__db =
            a_o__db + k_hat_1__db_per_meter * d__meter + k_hat_2__db_per_meter * d__meter.ln();
        *propmode = MODE_LINE_OF_SIGHT;
    } else {
        let d_5__meter = d_ML__meter + 200e3;
        let d_6__meter = d_ML__meter + 400e3;

        let mut h0 = -1.0f64;
        let a_6__db = troposcatter_loss(
            d_6__meter,
            theta_hzn,
            d_hzn__meter,
            h_e__meter,
            a_e__meter,
            n_s,
            f__mhz,
            theta_los,
            &mut h0,
        );
        let a_5__db = troposcatter_loss(
            d_5__meter,
            theta_hzn,
            d_hzn__meter,
            h_e__meter,
            a_e__meter,
            n_s,
            f__mhz,
            theta_los,
            &mut h0,
        );

        let m_s: f64;
        let a_s0__db: f64;
        let d_x__meter: f64;

        if a_5__db < 1000.0 {
            m_s = (a_6__db - a_5__db) / 200e3;

            d_x__meter = itm_max(
                itm_max(
                    d_sML__meter,
                    d_ML__meter
                        + 1.088 * (a_e__meter.powi(2) / f__mhz).powf(1.0 / 3.0) * f__mhz.ln(),
                ),
                (a_5__db - a_d0__db - m_s * d_5__meter) / (m_d - m_s),
            );

            a_s0__db = (m_d - m_s) * d_x__meter + a_d0__db;
        } else {
            m_s = m_d;
            a_s0__db = a_d0__db;
            d_x__meter = 10e6;
        }

        if d__meter > d_x__meter {
            *a_ref__db = m_s * d__meter + a_s0__db;
            *propmode = MODE_TROPOSCATTER;
        } else {
            *a_ref__db = m_d * d__meter + a_d0__db;
            *propmode = MODE_DIFFRACTION;
        }
    }

    *a_ref__db = itm_max(*a_ref__db, 0.0);

    SUCCESS
}

pub struct ItmOutput {
    pub a__db: f64,
    pub warnings: i32,
    pub inter_values: IntermediateValues,
}

pub fn itm_p2p_tls_ex(
    h_tx__meter: f64,
    h_rx__meter: f64,
    pfl: &[f64],
    climate: i32,
    n_0: f64,
    f__mhz: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    time: f64,
    location: f64,
    situation: f64,
) -> Result<ItmOutput, i32> {
    let mut inter_values = IntermediateValues::default();
    let mut warnings: i32 = NO_WARNINGS;

    let rtn = validate_inputs(
        h_tx__meter,
        h_rx__meter,
        climate,
        time,
        location,
        situation,
        n_0,
        f__mhz,
        pol,
        epsilon,
        sigma,
        mdvar,
        &mut warnings,
    );
    if rtn != SUCCESS {
        return Err(rtn);
    }

    inter_values.d__km = (pfl[0] * pfl[1]) / 1000.0;

    let np = pfl[0] as i32;

    let p10 = (0.1 * np as f64) as i32;
    let mut h_sys__meter = 0.0f64;

    for i in p10..=(np - p10) {
        h_sys__meter += pfl[i as usize + 2];
    }

    h_sys__meter /= (np - 2 * p10 + 1) as f64;

    let mut z_g = ComplexDouble::new(0.0, 0.0);
    let mut gamma_e = 0.0f64;
    let mut n_s = 0.0f64;

    initialize_point_to_point(
        f__mhz,
        h_sys__meter,
        n_0,
        pol,
        epsilon,
        sigma,
        &mut z_g,
        &mut gamma_e,
        &mut n_s,
    );

    let h__meter = [h_tx__meter, h_rx__meter];
    let mut theta_hzn = [0.0f64; 2];
    let mut d_hzn__meter = [0.0f64; 2];
    let mut h_e__meter = [0.0f64; 2];
    let mut delta_h__meter = 0.0f64;
    let mut d__meter = 0.0f64;

    quick_pfl(
        pfl,
        gamma_e,
        h__meter,
        &mut theta_hzn,
        &mut d_hzn__meter,
        &mut h_e__meter,
        &mut delta_h__meter,
        &mut d__meter,
    );

    let mut a_ref__db = 0.0f64;
    let mut propmode = MODE_NOT_SET;
    let rtn = longley_rice(
        theta_hzn,
        f__mhz,
        z_g,
        d_hzn__meter,
        h_e__meter,
        gamma_e,
        n_s,
        delta_h__meter,
        h__meter,
        d__meter,
        MODE_P2P,
        &mut a_ref__db,
        &mut warnings,
        &mut propmode,
    );
    if rtn != SUCCESS {
        return Err(rtn);
    }

    let a_fs__db = free_space_loss(d__meter, f__mhz);

    let a__db = variability(
        time,
        location,
        situation,
        h_e__meter,
        delta_h__meter,
        f__mhz,
        d__meter,
        a_ref__db,
        climate,
        mdvar,
        &mut warnings,
    ) + a_fs__db;

    inter_values.a_ref__db = a_ref__db;
    inter_values.a_fs__db = a_fs__db;
    inter_values.delta_h__meter = delta_h__meter;
    inter_values.d_hzn__meter = d_hzn__meter;
    inter_values.h_e__meter = h_e__meter;
    inter_values.n_s = n_s;
    inter_values.theta_hzn = theta_hzn;
    inter_values.mode = propmode;

    if warnings != NO_WARNINGS {
        Ok(ItmOutput {
            a__db,
            warnings: SUCCESS_WITH_WARNINGS,
            inter_values,
        })
    } else {
        Ok(ItmOutput {
            a__db,
            warnings: SUCCESS,
            inter_values,
        })
    }
}

pub fn itm_p2p_tls(
    h_tx__meter: f64,
    h_rx__meter: f64,
    pfl: &[f64],
    climate: i32,
    n_0: f64,
    f__mhz: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    time: f64,
    location: f64,
    situation: f64,
) -> Result<ItmOutput, i32> {
    itm_p2p_tls_ex(
        h_tx__meter,
        h_rx__meter,
        pfl,
        climate,
        n_0,
        f__mhz,
        pol,
        epsilon,
        sigma,
        mdvar,
        time,
        location,
        situation,
    )
}

pub fn itm_p2p_cr(
    h_tx__meter: f64,
    h_rx__meter: f64,
    pfl: &[f64],
    climate: i32,
    n_0: f64,
    f__mhz: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    confidence: f64,
    reliability: f64,
) -> Result<ItmOutput, i32> {
    let result = itm_p2p_tls_ex(
        h_tx__meter,
        h_rx__meter,
        pfl,
        climate,
        n_0,
        f__mhz,
        pol,
        epsilon,
        sigma,
        mdvar,
        reliability,
        50.0,
        confidence,
    );

    match result {
        Err(ERROR_INVALID_TIME) => Err(ERROR_INVALID_RELIABILITY),
        Err(ERROR_INVALID_SITUATION) => Err(ERROR_INVALID_CONFIDENCE),
        _ => result,
    }
}

pub fn itm_p2p_cr_ex(
    h_tx__meter: f64,
    h_rx__meter: f64,
    pfl: &[f64],
    climate: i32,
    n_0: f64,
    f__mhz: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    confidence: f64,
    reliability: f64,
) -> Result<ItmOutput, i32> {
    let result = itm_p2p_tls_ex(
        h_tx__meter,
        h_rx__meter,
        pfl,
        climate,
        n_0,
        f__mhz,
        pol,
        epsilon,
        sigma,
        mdvar,
        reliability,
        50.0,
        confidence,
    );

    match result {
        Err(ERROR_INVALID_TIME) => Err(ERROR_INVALID_RELIABILITY),
        Err(ERROR_INVALID_SITUATION) => Err(ERROR_INVALID_CONFIDENCE),
        _ => result,
    }
}

pub fn itm_area_tls_ex(
    h_tx__meter: f64,
    h_rx__meter: f64,
    tx_site_criteria: i32,
    rx_site_criteria: i32,
    d__km: f64,
    delta_h__meter: f64,
    climate: i32,
    n_0: f64,
    f__mhz: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    time: f64,
    location: f64,
    situation: f64,
) -> Result<ItmOutput, i32> {
    let mut inter_values = IntermediateValues::default();
    let mut warnings: i32 = NO_WARNINGS;

    let rtn = validate_inputs(
        h_tx__meter,
        h_rx__meter,
        climate,
        time,
        location,
        situation,
        n_0,
        f__mhz,
        pol,
        epsilon,
        sigma,
        mdvar,
        &mut warnings,
    );
    if rtn != SUCCESS {
        return Err(rtn);
    }

    if d__km <= 0.0 {
        return Err(ERROR_PATH_DISTANCE);
    }
    if delta_h__meter < 0.0 {
        return Err(ERROR_DELTA_H);
    }
    if tx_site_criteria != SITING_CRITERIA_RANDOM
        && tx_site_criteria != SITING_CRITERIA_CAREFUL
        && tx_site_criteria != SITING_CRITERIA_VERY_CAREFUL
    {
        return Err(ERROR_TX_SITING_CRITERIA);
    }
    if rx_site_criteria != SITING_CRITERIA_RANDOM
        && rx_site_criteria != SITING_CRITERIA_CAREFUL
        && rx_site_criteria != SITING_CRITERIA_VERY_CAREFUL
    {
        return Err(ERROR_RX_SITING_CRITERIA);
    }

    let site_criteria = [tx_site_criteria, rx_site_criteria];
    let h__meter = [h_tx__meter, h_rx__meter];
    inter_values.d__km = d__km;

    let mut theta_hzn = [0.0f64; 2];
    let mut d_hzn__meter = [0.0f64; 2];
    let mut h_e__meter = [0.0f64; 2];
    let mut z_g = ComplexDouble::new(0.0, 0.0);
    let mut n_s = 0.0f64;
    let mut gamma_e = 0.0f64;
    let mut a_ref__db = 0.0f64;

    initialize_point_to_point(
        f__mhz,
        0.0,
        n_0,
        pol,
        epsilon,
        sigma,
        &mut z_g,
        &mut gamma_e,
        &mut n_s,
    );

    initialize_area(
        site_criteria,
        gamma_e,
        delta_h__meter,
        h__meter,
        &mut h_e__meter,
        &mut d_hzn__meter,
        &mut theta_hzn,
    );

    let d__meter = d__km * 1000.0;
    let mut propmode = MODE_NOT_SET;
    let rtn = longley_rice(
        theta_hzn,
        f__mhz,
        z_g,
        d_hzn__meter,
        h_e__meter,
        gamma_e,
        n_s,
        delta_h__meter,
        h__meter,
        d__meter,
        MODE_AREA,
        &mut a_ref__db,
        &mut warnings,
        &mut propmode,
    );
    if rtn != SUCCESS {
        return Err(rtn);
    }

    let a_fs__db = free_space_loss(d__meter, f__mhz);

    let a__db = a_fs__db
        + variability(
            time,
            location,
            situation,
            h_e__meter,
            delta_h__meter,
            f__mhz,
            d__meter,
            a_ref__db,
            climate,
            mdvar,
            &mut warnings,
        );

    inter_values.a_ref__db = a_ref__db;
    inter_values.a_fs__db = a_fs__db;
    inter_values.delta_h__meter = delta_h__meter;
    inter_values.d_hzn__meter = d_hzn__meter;
    inter_values.h_e__meter = h_e__meter;
    inter_values.n_s = n_s;
    inter_values.theta_hzn = theta_hzn;
    inter_values.mode = propmode;

    if warnings != NO_WARNINGS {
        Ok(ItmOutput {
            a__db,
            warnings: SUCCESS_WITH_WARNINGS,
            inter_values,
        })
    } else {
        Ok(ItmOutput {
            a__db,
            warnings: SUCCESS,
            inter_values,
        })
    }
}

pub fn itm_area_tls(
    h_tx__meter: f64,
    h_rx__meter: f64,
    tx_site_criteria: i32,
    rx_site_criteria: i32,
    d__km: f64,
    delta_h__meter: f64,
    climate: i32,
    n_0: f64,
    f__mhz: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    time: f64,
    location: f64,
    situation: f64,
) -> Result<ItmOutput, i32> {
    itm_area_tls_ex(
        h_tx__meter,
        h_rx__meter,
        tx_site_criteria,
        rx_site_criteria,
        d__km,
        delta_h__meter,
        climate,
        n_0,
        f__mhz,
        pol,
        epsilon,
        sigma,
        mdvar,
        time,
        location,
        situation,
    )
}

pub fn itm_area_cr(
    h_tx__meter: f64,
    h_rx__meter: f64,
    tx_site_criteria: i32,
    rx_site_criteria: i32,
    d__km: f64,
    delta_h__meter: f64,
    climate: i32,
    n_0: f64,
    f__mhz: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    confidence: f64,
    reliability: f64,
) -> Result<ItmOutput, i32> {
    let result = itm_area_tls_ex(
        h_tx__meter,
        h_rx__meter,
        tx_site_criteria,
        rx_site_criteria,
        d__km,
        delta_h__meter,
        climate,
        n_0,
        f__mhz,
        pol,
        epsilon,
        sigma,
        mdvar,
        reliability,
        50.0,
        confidence,
    );

    match result {
        Err(ERROR_INVALID_TIME) => Err(ERROR_INVALID_RELIABILITY),
        Err(ERROR_INVALID_SITUATION) => Err(ERROR_INVALID_CONFIDENCE),
        _ => result,
    }
}

pub fn itm_area_cr_ex(
    h_tx__meter: f64,
    h_rx__meter: f64,
    tx_site_criteria: i32,
    rx_site_criteria: i32,
    d__km: f64,
    delta_h__meter: f64,
    climate: i32,
    n_0: f64,
    f__mhz: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    confidence: f64,
    reliability: f64,
) -> Result<ItmOutput, i32> {
    let result = itm_area_tls_ex(
        h_tx__meter,
        h_rx__meter,
        tx_site_criteria,
        rx_site_criteria,
        d__km,
        delta_h__meter,
        climate,
        n_0,
        f__mhz,
        pol,
        epsilon,
        sigma,
        mdvar,
        reliability,
        50.0,
        confidence,
    );

    match result {
        Err(ERROR_INVALID_TIME) => Err(ERROR_INVALID_RELIABILITY),
        Err(ERROR_INVALID_SITUATION) => Err(ERROR_INVALID_CONFIDENCE),
        _ => result,
    }
}
