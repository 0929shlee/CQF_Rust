use std::cmp::max;
use std::fs::File;
use std::io::Read;
use std::time::Instant;
use crate::model::CqiMatrix;
use crate::model::ConnectionMatrix;
use crate::model::comp_quality;

mod model;

fn main() {
    let info = read_info();
    let mut cqi_matrix = CqiMatrix::new(
        info.0,
        info.1,
        info.2,
        info.3,
        info.4,
        info.5,
        info.6,
    );
    cqi_matrix.generate();
    cqi_matrix.write();

    let mut connection_matrix = ConnectionMatrix::new(cqi_matrix);

    let now = Instant::now();

    match &info.7 as &str {
        "BRUTE_FORCE" => connection_matrix.brute_force(),
        "CQI_SORTING" => connection_matrix.cqi_sorting(),
        _ => panic!("Algorithm Name error"),
    }

    connection_matrix.write();

    println!("time spending: {}", now.elapsed().as_millis());
    let cq = comp_quality::calculate(&connection_matrix);
    println!("\n------------------------------");
    println!("CoMP Quality: {}", cq);
    println!("------------------------------");
    comp_quality::write(cq);

    /*
    test(100);
     */
}

fn read_info() -> (usize, usize, usize, usize, f64, f64, usize, String) {
    let mut file = match File::open("info.txt") {
        Err(e) => panic!("File open error: {}", e),
        Ok(file) => file,
    };

    let mut text = String::new();
    match file.read_to_string(&mut text) {
        Err(e) => panic!("file read to string error: {}", e),
        Ok(_) => {},
    }

    let str_arr = text.lines().collect::<Vec<&str>>();
    let n_gnb = str_arr[0]
        .split_whitespace()
        .collect::<Vec<&str>>()[1]
        .parse::<usize>()
        .unwrap();
    let n_ue = str_arr[1]
        .split_whitespace()
        .collect::<Vec<&str>>()[1]
        .parse::<usize>()
        .unwrap();
    let n_time = str_arr[2]
        .split_whitespace()
        .collect::<Vec<&str>>()[1]
        .parse::<usize>()
        .unwrap();
    let max_gnb_connection = str_arr[3]
        .split_whitespace()
        .collect::<Vec<&str>>()[1]
        .parse::<usize>()
        .unwrap();
    let scaling_factor = str_arr[4]
        .split_whitespace()
        .collect::<Vec<&str>>()[1]
        .parse::<f64>()
        .unwrap();
    let noise_density_factor = str_arr[5]
        .split_whitespace()
        .collect::<Vec<&str>>()[1]
        .parse::<f64>()
        .unwrap();
    let noise_length_factor = str_arr[6]
        .split_whitespace()
        .collect::<Vec<&str>>()[1]
        .parse::<usize>()
        .unwrap();
    let algorithm_name = str_arr[7]
        .split_whitespace()
        .collect::<Vec<&str>>()[1];


    (n_gnb,
     n_ue,
     n_time,
     max_gnb_connection,
     scaling_factor,
     noise_density_factor,
     noise_length_factor,
     String::from(algorithm_name)
    )
}

fn test(n_test: usize) {
    let mut error_cnt: usize = 0;
    for t in 0..n_test {
        println!("Running test... ({}/{})", t, n_test);
        let info = read_info();
        let mut cqi_matrix = CqiMatrix::new(
            info.0,
            info.1,
            info.2,
            info.3,
            info.4,
            info.5,
            info.6,
        );
        cqi_matrix.generate();

        let mut connection_matrix = ConnectionMatrix::new(cqi_matrix);

        let now = Instant::now();

        match &info.7 as &str {
            "BRUTE_FORCE" => connection_matrix.brute_force(),
            "CQI_SORTING" => connection_matrix.cqi_sorting(),
            _ => panic!("Algorithm Name error"),
        }

        if !connection_matrix.is_valid() {
            connection_matrix.write();
            assert!(false);
            error_cnt += 1;
        }

        println!("time spending: {}", now.elapsed().as_millis());
    }
    println!("Test failed for {} times", error_cnt);
}
