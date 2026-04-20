use crate::constants::{A_0_METER, MODE_P2P, THIRD};
use crate::helper::{fresnel_integral, height_function, itm_min};
use crate::terrain::{sigma_h_function, terrain_roughness};
use crate::types::ComplexDouble;

pub fn knife_edge_diffraction(
    d__meter: f64,
    f__mhz: f64,
    a_e__meter: f64,
    theta_los: f64,
    d_hzn__meter: [f64; 2],
) -> f64 {
    let d_ML__meter = d_hzn__meter[0] + d_hzn__meter[1];
    let theta_nlos = d__meter / a_e__meter - theta_los;

    let d_nlos__meter = d__meter - d_ML__meter;

    let v_1 = 0.0795775 * (f__mhz / 47.7) * theta_nlos.powi(2) * d_hzn__meter[0] * d_nlos__meter
        / (d_nlos__meter + d_hzn__meter[0]);
    let v_2 = 0.0795775 * (f__mhz / 47.7) * theta_nlos.powi(2) * d_hzn__meter[1] * d_nlos__meter
        / (d_nlos__meter + d_hzn__meter[1]);

    fresnel_integral(v_1) + fresnel_integral(v_2)
}

pub fn smooth_earth_diffraction(
    d__meter: f64,
    f__mhz: f64,
    a_e__meter: f64,
    theta_los: f64,
    d_hzn__meter: [f64; 2],
    h_e__meter: [f64; 2],
    z_g: ComplexDouble,
) -> f64 {
    let theta_nlos = d__meter / a_e__meter - theta_los;
    let d_ML__meter = d_hzn__meter[0] + d_hzn__meter[1];

    let mut a__meter = [0.0f64; 3];
    let mut d__km = [0.0f64; 3];
    let mut f_x__db = [0.0f64; 2];
    let mut k = [0.0f64; 3];
    let mut b_0 = [0.0f64; 3];
    let mut x__km = [0.0f64; 3];
    let mut c_0 = [0.0f64; 3];

    a__meter[0] = (d__meter - d_ML__meter) / (d__meter / a_e__meter - theta_los);
    a__meter[1] = 0.5 * d_hzn__meter[0].powi(2) / h_e__meter[0];
    a__meter[2] = 0.5 * d_hzn__meter[1].powi(2) / h_e__meter[1];

    d__km[0] = a__meter[0] * theta_nlos / 1000.0;
    d__km[1] = d_hzn__meter[0] / 1000.0;
    d__km[2] = d_hzn__meter[1] / 1000.0;

    for i in 0..3 {
        c_0[i] = ((4.0 / 3.0) * A_0_METER / a__meter[i]).powf(THIRD);
        k[i] = 0.017778 * c_0[i] * f__mhz.powf(-THIRD) / z_g.norm();
        b_0[i] = 1.607 - k[i];
    }

    x__km[1] = b_0[1] * c_0[1].powi(2) * f__mhz.powf(THIRD) * d__km[1];
    x__km[2] = b_0[2] * c_0[2].powi(2) * f__mhz.powf(THIRD) * d__km[2];
    x__km[0] = b_0[0] * c_0[0].powi(2) * f__mhz.powf(THIRD) * d__km[0] + x__km[1] + x__km[2];

    f_x__db[0] = height_function(x__km[1], k[1]);
    f_x__db[1] = height_function(x__km[2], k[2]);

    let g_x__db = 0.05751 * x__km[0] - 10.0 * x__km[0].log10();

    g_x__db - f_x__db[0] - f_x__db[1] - 20.0
}

pub fn diffraction_loss(
    d__meter: f64,
    d_hzn__meter: [f64; 2],
    h_e__meter: [f64; 2],
    z_g: ComplexDouble,
    a_e__meter: f64,
    delta_h__meter: f64,
    h__meter: [f64; 2],
    mode: i32,
    theta_los: f64,
    d_sML__meter: f64,
    f__mhz: f64,
) -> f64 {
    let a_k__db = knife_edge_diffraction(d__meter, f__mhz, a_e__meter, theta_los, d_hzn__meter);
    let a_se__db = smooth_earth_diffraction(
        d__meter,
        f__mhz,
        a_e__meter,
        theta_los,
        d_hzn__meter,
        h_e__meter,
        z_g,
    );

    let delta_h_dsML__meter = terrain_roughness(d_sML__meter, delta_h__meter);
    let sigma_h_d__meter = sigma_h_function(delta_h_dsML__meter);

    let a_fo__db = itm_min(
        15.0,
        5.0 * (1.0 + 1e-5 * h__meter[0] * h__meter[1] * f__mhz * sigma_h_d__meter).log10(),
    );

    let delta_h_d__meter = terrain_roughness(d__meter, delta_h__meter);

    let mut q = h__meter[0] * h__meter[1];
    let qk = h_e__meter[0] * h_e__meter[1] - q;

    if mode == MODE_P2P {
        q += 10.0;
    }

    let term1 = (1.0 + qk / q).sqrt();

    let d_ML__meter = d_hzn__meter[0] + d_hzn__meter[1];
    q = (term1 + (-theta_los * a_e__meter + d_ML__meter) / d__meter)
        * itm_min(delta_h_d__meter * f__mhz / 47.7, 6283.2);

    let w = 25.1 / (25.1 + q.sqrt());

    w * a_se__db + (1.0 - w) * a_k__db + a_fo__db
}
