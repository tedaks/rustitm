use crate::types::ComplexDouble;
use crate::constants::{GAMMA_A, PI};
use crate::enums::*;
use crate::helper::{itm_min, itm_max, fortran_dim};
use crate::terrain::{find_horizons, compute_delta_h, linear_least_squares_fit};

pub fn initialize_point_to_point(
    f__mhz: f64,
    h_sys__meter: f64,
    n_0: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    z_g: &mut ComplexDouble,
    gamma_e: &mut f64,
    n_s: &mut f64,
) {
    if h_sys__meter == 0.0 {
        *n_s = n_0;
    } else {
        *n_s = n_0 * (-h_sys__meter / 9460.0).exp();
    }

    *gamma_e = GAMMA_A * (1.0 - 0.04665 * (*n_s / 179.3).exp());

    let ep_r = ComplexDouble::new(epsilon, 18000.0 * sigma / f__mhz);

    *z_g = (ep_r - ComplexDouble::new(1.0, 0.0)).sqrt();

    if pol == POLARIZATION_VERTICAL {
        *z_g = *z_g / ep_r;
    }
}

pub fn initialize_area(
    site_criteria: [i32; 2],
    gamma_e: f64,
    delta_h__meter: f64,
    h__meter: [f64; 2],
    h_e__meter: &mut [f64; 2],
    d_hzn__meter: &mut [f64; 2],
    theta_hzn: &mut [f64; 2],
) {
    for i in 0..2 {
        if site_criteria[i] == SITING_CRITERIA_RANDOM {
            h_e__meter[i] = h__meter[i];
        } else {
            let mut b: f64;

            if site_criteria[i] == SITING_CRITERIA_CAREFUL {
                b = 4.0;
            } else {
                b = 9.0;
            }

            if h__meter[i] < 5.0 {
                b = b * (0.1 * PI * h__meter[i]).sin();
            }

            h_e__meter[i] = h__meter[i] + (1.0 + b) * (-itm_min(20.0, 2.0 * h__meter[i] / itm_max(1e-3, delta_h__meter))).exp();
        }

        let d_ls__meter = (2.0 * h_e__meter[i] / gamma_e).sqrt();

        let h_3__meter = 5.0;
        d_hzn__meter[i] = d_ls__meter * (-0.07 * (delta_h__meter / itm_max(h_e__meter[i], h_3__meter)).sqrt()).exp();

        theta_hzn[i] = (0.65 * delta_h__meter * (d_ls__meter / d_hzn__meter[i] - 1.0) - 2.0 * h_e__meter[i]) / d_ls__meter;
    }
}

pub fn quick_pfl(
    pfl: &[f64],
    gamma_e: f64,
    h__meter: [f64; 2],
    theta_hzn: &mut [f64; 2],
    d_hzn__meter: &mut [f64; 2],
    h_e__meter: &mut [f64; 2],
    delta_h__meter: &mut f64,
    d__meter: &mut f64,
) {
    let mut fit_tx: f64 = 0.0;
    let mut fit_rx: f64 = 0.0;
    let mut q: f64;
    let d_start__meter: f64;
    let d_end__meter: f64;

    *d__meter = pfl[0] * pfl[1];

    let np = pfl[0] as i32;

    let a_e__meter = 1.0 / gamma_e;

    find_horizons(pfl, a_e__meter, h__meter, theta_hzn, d_hzn__meter);

    d_start__meter = itm_min(15.0 * h__meter[0], 0.1 * d_hzn__meter[0]);
    d_end__meter = *d__meter - itm_min(15.0 * h__meter[1], 0.1 * d_hzn__meter[1]);

    *delta_h__meter = compute_delta_h(pfl, d_start__meter, d_end__meter);

    if d_hzn__meter[0] + d_hzn__meter[1] > 1.5 * *d__meter {
        linear_least_squares_fit(pfl, d_start__meter, d_end__meter, &mut fit_tx, &mut fit_rx);

        h_e__meter[0] = h__meter[0] + fortran_dim(pfl[2], fit_tx);
        h_e__meter[1] = h__meter[1] + fortran_dim(pfl[(np + 2) as usize], fit_rx);

        for i in 0..2 {
            d_hzn__meter[i] = (2.0 * h_e__meter[i] * a_e__meter).sqrt() * (-0.07 * (*delta_h__meter / itm_max(h_e__meter[i], 5.0)).sqrt()).exp();
        }

        let combined_horizons__meter = d_hzn__meter[0] + d_hzn__meter[1];
        if combined_horizons__meter <= *d__meter {
            q = (*d__meter / combined_horizons__meter).powi(2);

            for i in 0..2 {
                h_e__meter[i] = h_e__meter[i] * q;
                d_hzn__meter[i] = (2.0 * h_e__meter[i] * a_e__meter).sqrt() * (-0.07 * (*delta_h__meter / itm_max(h_e__meter[i], 5.0)).sqrt()).exp();
            }
        }

        for i in 0..2 {
            q = (2.0 * h_e__meter[i] * a_e__meter).sqrt();
            theta_hzn[i] = (0.65 * *delta_h__meter * (q / d_hzn__meter[i] - 1.0) - 2.0 * h_e__meter[i]) / q;
        }
    } else {
        let mut dummy = 0.0;

        linear_least_squares_fit(pfl, d_start__meter, 0.9 * d_hzn__meter[0], &mut fit_tx, &mut dummy);
        h_e__meter[0] = h__meter[0] + fortran_dim(pfl[2], fit_tx);

        linear_least_squares_fit(pfl, *d__meter - 0.9 * d_hzn__meter[1], d_end__meter, &mut dummy, &mut fit_rx);
        h_e__meter[1] = h__meter[1] + fortran_dim(pfl[(np + 2) as usize], fit_rx);
    }
}