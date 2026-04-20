use rustitm::entry::itm_p2p_tls;

fn main() {
    let pfl = vec![
        10.0, 100.0, 100.0, 110.0, 120.0, 115.0, 105.0, 100.0, 95.0, 90.0, 85.0, 80.0,
    ];

    match itm_p2p_tls(
        100.0, 10.0, &pfl, 5, 301.0, 1000.0, 1, 15.0, 0.001, 3, 50.0, 50.0, 50.0,
    ) {
        Ok(output) => println!(
            "Attenuation: {:.2} dB, Warnings: {}",
            output.a__db, output.warnings
        ),
        Err(code) => println!("Error: {}", code),
    }
}
