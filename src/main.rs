use std::fs::File;
use std::io::{Read, stdout, Write};
use std::path::Path;
use std::time::Instant;
use crate::model::CqiMatrix;
use crate::model::ConnectionMatrix;
use crate::model::comp_quality;

mod model;

type InfoType = (usize, usize, usize, usize, f64, f64, usize, f64, String);

fn main() {
    let info = read_info();
    /*
    match &info.8 as &str {
        "BF" => run_brute_force(&info),
        "CS" => run_cqi_sorting(&info),
        "ES" => run_expectation_sorting(&info),
        "SS" => run_shuffle_split(&info, 100),
        _ => panic!("Algorithm Name error"),
    }
     */

    for i in 2..11 {
        /*
        let info = (info.0, info.1, info.2, i as usize, info.4, info.5, info.6, info.7, String::from(&info.8));
        run_cqi_sorting(&info);
        run_expectation_sorting(&info);
        run_naive_approach(&info);

         */

        /*
        let info = (info.0, info.1, info.2, i * 10 as usize, info.4, info.5, info.6, info.7, String::from(&info.8));
        run_cqi_sorting(&info);
        run_expectation_sorting(&info);
        run_naive_approach(&info);
         */
    }
    for i in 2..11 {
        let info = (info.0, info.1, info.2, i as usize, info.4, info.5, info.6, info.7, String::from(&info.8));
        run_brute_force(&info);
    }

    /*
     */

    /*
    test(100);
     */
}


fn run_expectation_sorting(info: &InfoType) {
    let mut cqi_matrix = get_cqi_matrix(info);
    cqi_matrix.generate();
    let mut connection_matrix = ConnectionMatrix::new(cqi_matrix);
    let now = Instant::now();
    connection_matrix.expectation_sorting();
    let rt = now.elapsed().as_millis();
    let cq = comp_quality::calculate(&connection_matrix);
    connection_matrix.write();
    write_results(format!("result_{}_{}.txt", connection_matrix.algorithm_name, info.3).as_str(),
                  cq, rt);
}
fn run_naive_approach(info: &InfoType) {
    let mut cqi_matrix = get_cqi_matrix(info);
    cqi_matrix.generate();
    let mut connection_matrix = ConnectionMatrix::new(cqi_matrix);
    let now = Instant::now();
    connection_matrix.naive_approach();
    let rt = now.elapsed().as_millis();
    let cq = comp_quality::calculate(&connection_matrix);
    connection_matrix.write();
    write_results(format!("result_{}_{}.txt", connection_matrix.algorithm_name, info.3).as_str(),
                  cq, rt);

}
fn run_brute_force(info: &InfoType) {
    let mut cqi_matrix = get_cqi_matrix(info);
    cqi_matrix.generate();
    let mut connection_matrix = ConnectionMatrix::new(cqi_matrix);
    let now = Instant::now();
    connection_matrix.brute_force();
    let rt = now.elapsed().as_millis();
    let cq = comp_quality::calculate(&connection_matrix);
    connection_matrix.write();
    write_results(format!("result_{}_{}.txt", connection_matrix.algorithm_name, info.3).as_str(),
                  cq, rt);
}
fn run_cqi_sorting(info: &InfoType) {
    let mut cqi_matrix = get_cqi_matrix(info);
    cqi_matrix.generate();
    let mut connection_matrix = ConnectionMatrix::new(cqi_matrix);
    let now = Instant::now();
    connection_matrix.cqi_sorting();
    let rt = now.elapsed().as_millis();
    let cq = comp_quality::calculate(&connection_matrix);
    connection_matrix.write();
    write_results(format!("result_{}_{}.txt", connection_matrix.algorithm_name, info.3).as_str(),
                  cq, rt);
}
fn run_shuffle_split(info: &InfoType, repeat: u64) {
    let mut cqi_matrix = get_cqi_matrix(info);
    cqi_matrix.generate();
    let mut connection_matrix = ConnectionMatrix::new(cqi_matrix);
    let now = Instant::now();
    connection_matrix.shuffle_split(repeat);
    let rt = now.elapsed().as_millis();
    let cq = comp_quality::calculate(&connection_matrix);
    connection_matrix.write();
    write_results(format!("result_{}_{}_{}.txt",
                          connection_matrix.algorithm_name,
                          repeat,
                          info.3).as_str(),
                  cq, rt);
}
fn get_cqi_matrix(info: &InfoType) -> CqiMatrix {
    CqiMatrix::new(
        info.0,
        info.1,
        info.2,
        info.3,
        info.4,
        info.5,
        info.6,
        info.7,
    )
}

fn write_results(file_path: &str, comp_quality: f64, running_time: u128) {
    let path = Path::new(file_path);
    let mut file = match File::create(&path) {
        Err(e) => panic!("File create error: {}", e),
        Ok(file) => file,
    };

    let s = format!("COMP_QUALITY: {}\nRUNNING_TIME: {}\n", comp_quality, running_time);
    match file.write_all(s.as_bytes()) {
        Err(e) => panic!("file write all error: {}", e),
        Ok(_) => {},
    }
}

fn read_info() -> InfoType {
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
    let penalty_factor = str_arr[7]
        .split_whitespace()
        .collect::<Vec<&str>>()[1]
        .parse::<f64>()
        .unwrap();
    let algorithm_name = str_arr[8]
        .split_whitespace()
        .collect::<Vec<&str>>()[1];

    (n_gnb,
     n_ue,
     n_time,
     max_gnb_connection,
     scaling_factor,
     noise_density_factor,
     noise_length_factor,
     penalty_factor,
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
            info.7,
        );
        cqi_matrix.generate();

        let mut connection_matrix = ConnectionMatrix::new(cqi_matrix);

        let now = Instant::now();

        match &info.8 as &str {
            "BRUTE_FORCE" => connection_matrix.brute_force(),
            "CQI_SORTING" => connection_matrix.cqi_sorting(),
            "SHUFFLE_SPLIT" => connection_matrix.shuffle_split(100),
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
