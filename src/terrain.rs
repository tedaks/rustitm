use crate::helper::fortran_dim;

pub fn terrain_roughness(d__meter: f64, delta_h__meter: f64) -> f64 {
    delta_h__meter * (1.0 - 0.8 * (-d__meter / 50e3).exp())
}

pub fn sigma_h_function(delta_h__meter: f64) -> f64 {
    0.78 * delta_h__meter * (-0.5 * delta_h__meter.powf(0.25)).exp()
}

pub fn compute_delta_h(pfl: &[f64], d_start__meter: f64, d_end__meter: f64) -> f64 {
    let mut s = [0.0f64; 247];

    let np = pfl[0] as i32;
    let xi = pfl[1];
    let mut x_start = d_start__meter / xi;
    let mut x_end = d_end__meter / xi;

    if x_end - x_start < 2.0 {
        return 0.0;
    }

    let mut p10 = (0.1 * (x_end - x_start + 8.0)) as i32;
    p10 = p10.max(4).min(25);

    let n = 10 * p10 - 5;
    let p90 = n - p10;

    let np_s = (n - 1) as f64;
    s[0] = np_s;
    s[1] = 1.0;

    x_end = (x_end - x_start) / np_s;
    let mut i = x_start as i32;
    x_start -= (i + 1) as f64;

    for j in 0..n {
        while x_start > 0.0 && (i + 1) < np {
            x_start -= 1.0;
            i += 1;
        }

        s[j as usize + 2] =
            pfl[i as usize + 3] + (pfl[i as usize + 3] - pfl[i as usize + 2]) * x_start;

        x_start += x_end;
    }

    let mut fit_y1 = 0.0f64;
    let mut fit_y2 = 0.0f64;
    linear_least_squares_fit(&s, 0.0, np_s, &mut fit_y1, &mut fit_y2);

    fit_y2 = (fit_y2 - fit_y1) / np_s;

    let mut diffs: Vec<f64> = Vec::with_capacity(n as usize);
    for j in 0..n {
        diffs.push(s[j as usize + 2] - fit_y1);
        fit_y1 += fit_y2;
    }

    diffs.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let q10 = diffs[p10 as usize - 1];
    let q90 = diffs[p90 as usize];

    let delta_h_d__meter = q10 - q90;

    let delta_h__meter =
        delta_h_d__meter / (1.0 - 0.8 * (-(d_end__meter - d_start__meter) / 50e3).exp());

    delta_h__meter
}

pub fn find_horizons(
    pfl: &[f64],
    a_e__meter: f64,
    h__meter: [f64; 2],
    theta_hzn: &mut [f64; 2],
    d_hzn__meter: &mut [f64; 2],
) {
    let np = pfl[0] as i32;
    let xi = pfl[1];
    let d__meter = pfl[0] * pfl[1];

    let z_tx__meter = pfl[2] + h__meter[0];
    let z_rx__meter = pfl[(np + 2) as usize] + h__meter[1];

    theta_hzn[0] = (z_rx__meter - z_tx__meter) / d__meter - d__meter / (2.0 * a_e__meter);
    theta_hzn[1] = -(z_rx__meter - z_tx__meter) / d__meter - d__meter / (2.0 * a_e__meter);

    d_hzn__meter[0] = d__meter;
    d_hzn__meter[1] = d__meter;

    let mut d_tx__meter = 0.0f64;
    let mut d_rx__meter = d__meter;

    for i in 1..np {
        d_tx__meter += xi;
        d_rx__meter -= xi;

        let theta_tx =
            (pfl[(i + 2) as usize] - z_tx__meter) / d_tx__meter - d_tx__meter / (2.0 * a_e__meter);
        let theta_rx =
            -(z_rx__meter - pfl[(i + 2) as usize]) / d_rx__meter - d_rx__meter / (2.0 * a_e__meter);

        if theta_tx > theta_hzn[0] {
            theta_hzn[0] = theta_tx;
            d_hzn__meter[0] = d_tx__meter;
        }

        if theta_rx > theta_hzn[1] {
            theta_hzn[1] = theta_rx;
            d_hzn__meter[1] = d_rx__meter;
        }
    }
}

pub fn linear_least_squares_fit(
    pfl: &[f64],
    d_start: f64,
    d_end: f64,
    fit_y1: &mut f64,
    fit_y2: &mut f64,
) {
    let np = pfl[0] as i32;

    let mut i_start = fortran_dim(d_start / pfl[1], 0.0) as i32;
    let mut i_end = np - fortran_dim(np as f64, d_end / pfl[1]) as i32;

    if i_end <= i_start {
        i_start = fortran_dim(i_start as f64, 1.0) as i32;
        i_end = np - fortran_dim(np as f64, i_end as f64 + 1.0) as i32;
    }

    let x_length = (i_end - i_start) as f64;

    let mut mid_shifted_index = -0.5 * x_length;
    let mid_shifted_end = i_end as f64 + mid_shifted_index;

    let mut sum_y = 0.5 * (pfl[i_start as usize + 2] + pfl[i_end as usize + 2]);
    let mut scaled_sum_y =
        0.5 * (pfl[i_start as usize + 2] - pfl[i_end as usize + 2]) * mid_shifted_index;

    let mut i_start_iter = i_start;
    for _i in 2..=(x_length as i32) {
        i_start_iter += 1;
        mid_shifted_index += 1.0;

        sum_y += pfl[i_start_iter as usize + 2];
        scaled_sum_y += pfl[i_start_iter as usize + 2] * mid_shifted_index;
    }

    sum_y = sum_y / x_length;
    scaled_sum_y = scaled_sum_y * 12.0 / ((x_length * x_length + 2.0) * x_length);

    *fit_y1 = sum_y - scaled_sum_y * mid_shifted_end;
    *fit_y2 = sum_y + scaled_sum_y * (np as f64 - mid_shifted_end);
}
