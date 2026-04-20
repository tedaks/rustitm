use crate::constants::{A_9000_METER, THIRD};
use crate::errors::*;
use crate::helper::{curve, inverse_complementary_cumulative_distribution_function};
use crate::terrain::terrain_roughness;

pub fn variability(
    time: f64,
    location: f64,
    situation: f64,
    h_e__meter: [f64; 2],
    delta_h__meter: f64,
    f__mhz: f64,
    d__meter: f64,
    a_ref__db: f64,
    climate: i32,
    mdvar: i32,
    warnings: &mut i32,
) -> f64 {
    let all_year: [[f64; 7]; 5] = [
        [-9.67, -0.62, 1.26, -9.21, -0.62, -0.39, 3.15],
        [12.7, 9.19, 15.5, 9.05, 9.19, 2.86, 857.9],
        [
            144.9e3, 228.9e3, 262.6e3, 84.1e3, 228.9e3, 141.7e3, 2222.0e3,
        ],
        [
            190.3e3, 205.2e3, 185.2e3, 101.1e3, 205.2e3, 315.9e3, 164.8e3,
        ],
        [133.8e3, 143.6e3, 99.8e3, 98.6e3, 143.6e3, 167.4e3, 116.3e3],
    ];

    let bsm1 = [2.13, 2.66, 6.11, 1.98, 2.68, 6.86, 8.51];
    let bsm2 = [159.5, 7.67, 6.65, 13.11, 7.16, 10.38, 169.8];
    let xsm1 = [762.2e3, 100.4e3, 138.2e3, 139.1e3, 93.7e3, 187.8e3, 609.8e3];
    let xsm2 = [
        123.6e3, 172.5e3, 242.2e3, 132.7e3, 186.8e3, 169.6e3, 119.9e3,
    ];
    let xsm3 = [94.5e3, 136.4e3, 178.6e3, 193.5e3, 133.5e3, 108.9e3, 106.6e3];

    let bsp1 = [2.11, 6.87, 10.08, 3.68, 4.75, 8.58, 8.43];
    let bsp2 = [102.3, 15.53, 9.60, 159.3, 8.12, 13.97, 8.19];
    let xsp1 = [636.9e3, 138.7e3, 165.3e3, 464.4e3, 93.2e3, 216.0e3, 136.2e3];
    let xsp2 = [134.8e3, 143.7e3, 225.7e3, 93.1e3, 135.9e3, 152.0e3, 188.5e3];
    let xsp3 = [95.6e3, 98.6e3, 129.7e3, 94.2e3, 113.4e3, 122.7e3, 122.9e3];

    let c_d = [1.224, 0.801, 1.380, 1.000, 1.224, 1.518, 1.518];
    let z_d = [1.282, 2.161, 1.282, 20.0, 1.282, 1.282, 1.282];

    let bfm1 = [1.0, 1.0, 1.0, 1.0, 0.92, 1.0, 1.0];
    let bfm2 = [0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0];
    let bfm3 = [0.0, 0.0, 0.0, 0.0, 1.77, 0.0, 0.0];

    let bfp1 = [1.0, 0.93, 1.0, 0.93, 0.93, 1.0, 1.0];
    let bfp2 = [0.0, 0.31, 0.0, 0.19, 0.31, 0.0, 0.0];
    let bfp3 = [0.0, 2.00, 0.0, 1.79, 2.00, 0.0, 0.0];

    let mut z_t = inverse_complementary_cumulative_distribution_function(time / 100.0);
    let mut z_l = inverse_complementary_cumulative_distribution_function(location / 100.0);
    let z_s = inverse_complementary_cumulative_distribution_function(situation / 100.0);

    let climate_idx = (climate - 1) as usize;

    let wn = f__mhz / 47.7;

    let d_ex__meter = (2.0 * A_9000_METER * h_e__meter[0]).sqrt()
        + (2.0 * A_9000_METER * h_e__meter[1]).sqrt()
        + (575.7e12 / wn).powf(THIRD);

    let d_e__meter;
    if d__meter < d_ex__meter {
        d_e__meter = 130e3 * d__meter / d_ex__meter;
    } else {
        d_e__meter = 130e3 + d__meter - d_ex__meter;
    }

    let mut mdvar_internal = mdvar;
    let plus20 = mdvar_internal >= 20;
    if plus20 {
        mdvar_internal -= 20;
    }

    let sigma_s;
    if plus20 {
        sigma_s = 0.0;
    } else {
        let d__meter = 100e3;
        sigma_s = 5.0 + 3.0 * (-d_e__meter / d__meter).exp();
    }

    let plus10 = mdvar_internal >= 10;
    if plus10 {
        mdvar_internal -= 10;
    }

    let v_med__db = curve(
        all_year[0][climate_idx],
        all_year[1][climate_idx],
        all_year[2][climate_idx],
        all_year[3][climate_idx],
        all_year[4][climate_idx],
        d_e__meter,
    );

    if mdvar_internal == 0 {
        z_t = z_s;
        z_l = z_s;
    } else if mdvar_internal == 1 {
        z_l = z_s;
    } else if mdvar_internal == 2 {
        z_l = z_t;
    }

    if z_t.abs() > 3.10 || z_l.abs() > 3.10 || z_s.abs() > 3.10 {
        *warnings |= WARN_EXTREME_VARIABILITIES;
    }

    let sigma_l;
    if plus10 {
        sigma_l = 0.0;
    } else {
        let delta_h_d__meter = terrain_roughness(d__meter, delta_h__meter);
        sigma_l = 10.0 * wn * delta_h_d__meter / (wn * delta_h_d__meter + 13.0);
    }
    let y_l = sigma_l * z_l;

    let q = (0.133 * wn).ln();
    let g_minus = bfm1[climate_idx] + bfm2[climate_idx] / ((bfm3[climate_idx] * q).powi(2) + 1.0);
    let g_plus = bfp1[climate_idx] + bfp2[climate_idx] / ((bfp3[climate_idx] * q).powi(2) + 1.0);

    let sigma_t_minus = curve(
        bsm1[climate_idx],
        bsm2[climate_idx],
        xsm1[climate_idx],
        xsm2[climate_idx],
        xsm3[climate_idx],
        d_e__meter,
    ) * g_minus;
    let sigma_t_plus = curve(
        bsp1[climate_idx],
        bsp2[climate_idx],
        xsp1[climate_idx],
        xsp2[climate_idx],
        xsp3[climate_idx],
        d_e__meter,
    ) * g_plus;

    let sigma_td = c_d[climate_idx] * sigma_t_plus;
    let tgtd = (sigma_t_plus - sigma_td) * z_d[climate_idx];

    let sigma_t;
    if z_t < 0.0 {
        sigma_t = sigma_t_minus;
    } else if z_t <= z_d[climate_idx] {
        sigma_t = sigma_t_plus;
    } else {
        sigma_t = sigma_td + tgtd / z_t;
    }
    let y_t = sigma_t * z_t;

    let y_s_temp =
        sigma_s.powi(2) + y_t.powi(2) / (7.8 + z_s.powi(2)) + y_l.powi(2) / (24.0 + z_s.powi(2));

    let (y_r, y_s): (f64, f64) = if mdvar_internal == 0 {
        (
            0.0,
            (sigma_t.powi(2) + sigma_l.powi(2) + y_s_temp).sqrt() * z_s,
        )
    } else if mdvar_internal == 1 {
        (y_t, (sigma_l.powi(2) + y_s_temp).sqrt() * z_s)
    } else if mdvar_internal == 2 {
        (
            (sigma_t.powi(2) + sigma_l.powi(2)).sqrt() * z_t,
            y_s_temp.sqrt() * z_s,
        )
    } else {
        (y_t + y_l, y_s_temp.sqrt() * z_s)
    };

    let mut result = a_ref__db - v_med__db - y_r - y_s;

    if result < 0.0 {
        result = result * (29.0 - result) / (29.0 - 10.0 * result);
    }

    result
}

pub fn free_space_loss(d__meter: f64, f__mhz: f64) -> f64 {
    32.45 + 20.0 * f__mhz.log10() + 20.0 * (d__meter / 1000.0).log10()
}

pub fn validate_inputs(
    h_tx__meter: f64,
    h_rx__meter: f64,
    climate: i32,
    time: f64,
    location: f64,
    situation: f64,
    n_0: f64,
    f__mhz: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    warnings: &mut i32,
) -> i32 {
    if h_tx__meter < 1.0 || h_tx__meter > 1000.0 {
        *warnings |= WARN_TX_TERMINAL_HEIGHT;
    }
    if h_tx__meter < 0.5 || h_tx__meter > 3000.0 {
        return ERROR_TX_TERMINAL_HEIGHT;
    }

    if h_rx__meter < 1.0 || h_rx__meter > 1000.0 {
        *warnings |= WARN_RX_TERMINAL_HEIGHT;
    }
    if h_rx__meter < 0.5 || h_rx__meter > 3000.0 {
        return ERROR_RX_TERMINAL_HEIGHT;
    }

    if climate != 1
        && climate != 2
        && climate != 3
        && climate != 4
        && climate != 5
        && climate != 6
        && climate != 7
    {
        return ERROR_INVALID_RADIO_CLIMATE;
    }

    if n_0 < 250.0 || n_0 > 400.0 {
        return ERROR_REFRACTIVITY;
    }

    if f__mhz < 40.0 || f__mhz > 10000.0 {
        *warnings |= WARN_FREQUENCY;
    }
    if f__mhz < 20.0 || f__mhz > 20000.0 {
        return ERROR_FREQUENCY;
    }

    if pol != 0 && pol != 1 {
        return ERROR_POLARIZATION;
    }

    if epsilon < 1.0 {
        return ERROR_EPSILON;
    }

    if sigma <= 0.0 {
        return ERROR_SIGMA;
    }

    if (mdvar < 0)
        || (mdvar > 3 && mdvar < 10)
        || (mdvar > 13 && mdvar < 20)
        || (mdvar > 23 && mdvar < 30)
        || (mdvar > 33)
    {
        return ERROR_MDVAR;
    }

    if situation <= 0.0 || situation >= 100.0 {
        return ERROR_INVALID_SITUATION;
    }

    if time <= 0.0 || time >= 100.0 {
        return ERROR_INVALID_TIME;
    }

    if location <= 0.0 || location >= 100.0 {
        return ERROR_INVALID_LOCATION;
    }

    SUCCESS
}
