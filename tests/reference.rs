//! Validation against the NTIA/ITS reference vectors in ../itm/.
//!
//! - p2p: cmd_examples (i_p2ptls.txt + pfl.txt -> o_p2ptls.txt, 114.5 dB)
//! - area: all 5 rows of area.csv

use rustitm::entry::{itm_area_tls, itm_p2p_tls};

const ITM_DIR: &str = "/home/bortre/02-lab/sources/itm";

fn read_csv_row(line: &str) -> Vec<f64> {
    line.trim()
        .split(',')
        .map(|s| s.trim().parse::<f64>().unwrap())
        .collect()
}

fn load_pfl(path: &str) -> Vec<f64> {
    let raw = std::fs::read_to_string(path).expect("read pfl");
    read_csv_row(&raw)
}

#[test]
fn p2p_tls_cmd_example_matches_ref() {
    let pfl = load_pfl(&format!("{}/cmd_examples/pfl.txt", ITM_DIR));

    // i_p2ptls.txt
    let out = itm_p2p_tls(
        15.0, 3.0, &pfl,
        5,      // climate: continental temperate
        301.0,  // N_0
        3500.0, // f__mhz
        1,      // pol: vertical
        15.0, 0.005,
        1,      // mdvar: accidental
        50.0, 50.0, 50.0,
    ).expect("itm_p2p_tls ok");

    // o_p2ptls.txt -> Basic Transmission Loss 114.5 dB
    let expected = 114.5;
    assert!(
        (out.a__db - expected).abs() < 0.1,
        "p2p cmd_example: got {:.4} dB, expected {:.4} dB",
        out.a__db,
        expected,
    );
}

#[test]
fn area_tls_csv_matches_ref() {
    let csv = std::fs::read_to_string(format!("{}/area.csv", ITM_DIR)).unwrap();
    let mut lines = csv.lines();
    lines.next(); // header

    let mut failures: Vec<String> = Vec::new();

    for (idx, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let r = read_csv_row(line);
        // h_tx,h_rx,delta_h,mdvar,d__km,tx_site,rx_site,epsilon,sigma,N_0,f_mhz,pol,climate,time,location,situation,A__db
        let (h_tx, h_rx, delta_h, mdvar, d_km, tx_site, rx_site,
             epsilon, sigma, n_0, f_mhz, pol, climate,
             time, location, situation, expected)
             = (r[0], r[1], r[2], r[3] as i32, r[4],
                r[5] as i32, r[6] as i32, r[7], r[8], r[9], r[10],
                r[11] as i32, r[12] as i32, r[13], r[14], r[15], r[16]);

        let got = itm_area_tls(
            h_tx, h_rx, tx_site, rx_site, d_km, delta_h,
            climate, n_0, f_mhz, pol, epsilon, sigma, mdvar,
            time, location, situation,
        );

        match got {
            Ok(out) => {
                let diff = (out.a__db - expected).abs();
                if diff >= 0.1 {
                    failures.push(format!(
                        "row {}: got {:.4} dB, expected {:.4} dB (diff {:.4})",
                        idx + 2, out.a__db, expected, diff,
                    ));
                }
            }
            Err(code) => {
                failures.push(format!("row {}: error code {}", idx + 2, code));
            }
        }
    }

    assert!(
        failures.is_empty(),
        "area.csv mismatches:\n  {}",
        failures.join("\n  "),
    );
}
