use std::fs::File;
use std::io::{Read, stdout, Write};
use std::path::Path;

pub struct Matrix {
    pub n_gnb: usize,
    pub n_ue: usize,
    pub n_time: usize,

    pub file_path: String,

    pub max_cqi: u8,
    pub max_gnb_connection: usize,

    pub matrix: Vec<Vec<Vec<u8>>>,
}

impl Matrix {
    pub fn new(
        n_gnb: usize,
        n_ue: usize,
        n_time: usize,
        file_path: String,
        max_gnb_connection: usize) -> Self {
        Matrix {
            n_gnb,
            n_ue,
            n_time,
            file_path,
            max_cqi: 15,
            max_gnb_connection,
            matrix: vec![vec![vec![0; n_time]; n_ue]; n_gnb],
        }
    }

    pub fn read(&mut self) {
        let mut file = match File::open(&self.file_path) {
            Err(e) => panic!("File open error: {}", e),
            Ok(file) => file,
        };

        self.file_to_matrix(&mut file);
    }
    pub fn write(&self) {
        let path = Path::new(&self.file_path);
        let mut file = match File::create(&path) {
            Err(e) => panic!("File create error: {}", e),
            Ok(file) => file,
        };

        self.matrix_to_file(&mut file);
    }
    pub fn is_there_file(&self) -> bool {
        match File::open(&self.file_path) {
            Err(_) => false,
            Ok(_) => true,
        }
    }
}

impl Matrix {
    //File read
    fn file_to_matrix(&mut self, file: &mut File) {
        let mut s = String::new();
        match file.read_to_string(&mut s) {
            Err(e) => panic!("file read to string error: {}", e),
            Ok(_) => {},
        }

        self.string_to_matrix(&s);
    }
    fn string_to_matrix(&mut self, s: &String) {
        self.update_info_from_string(&s);

        let mut cnt_time: usize = 0;
        let mut cnt_gnb: usize = 0;
        for s in s.lines().collect::<Vec<&str>>() {
            if s != "" {
                let mut cnt_ue: usize = 0;
                let str_num_arr = s.split_whitespace().collect::<Vec<&str>>();
                for str_num in str_num_arr {
                    self.matrix[cnt_gnb][cnt_ue][cnt_time] = str_num.parse::<u8>().unwrap();
                    cnt_ue += 1;
                }
                cnt_gnb += 1;
            }
            else {
                cnt_time += 1;
                cnt_gnb = 0;
            }
        }
    }
    fn update_info_from_string(&mut self, s: &String) {
        let str_arr = s.lines().collect::<Vec<&str>>();

        let n_time = str_arr.iter().filter(|str| **str == "").count();
        let n_gnb = str_arr.iter().count() / n_time - 1;
        let n_ue = str_arr[0].split_whitespace().count();

        self.n_gnb = n_gnb;
        self.n_ue = n_ue;
        self.n_time = n_time;
        self.matrix = vec![vec![vec![0; n_time]; n_ue]; n_gnb];
    }

}

impl Matrix {
    //File write
    fn matrix_to_file(&self, file: &mut File) {
        let mut s = String::new();
        self.matrix_to_string(&mut s);
        match file.write_all(s.as_bytes()) {
            Err(e) => panic!("file write all error: {}", e),
            Ok(_) => {},
        }
    }
    fn matrix_to_string(&self, s: &mut String) {
        for t in 0..self.n_time {
            self.matrix_of_time_to_string(t, s);
            s.push_str("\n");
        }
    }
    fn matrix_of_time_to_string(&self, time_idx: usize, s: &mut String) {
        for g in 0..self.n_gnb {
            for u in 0..self.n_ue {
                let str_num = format!("{0: <3}", self.matrix[g][u][time_idx]);
                s.push_str(&str_num);
            }
            s.push_str("\n");
        }
    }
}

impl Matrix {
    //Print
    pub fn print(&self) {
        for t in 0..self.n_time {
            self._print_of_time(t);
            println!();
        }
    }
    pub fn print_of_time(&self, time_idx: usize) {
        println!("----------TIME: {}----------", time_idx + 1);
        self._print_of_time(time_idx);
        println!("----------------------------");
    }
    fn _print_of_time(&self, time_idx: usize) {
        for g in 0..self.n_gnb {
            for u in 0..self.n_ue {
                if self.matrix[g][u][time_idx] == 0 {
                    print!("-  ");
                }
                else {
                    print!("{0: <3}", self.matrix[g][u][time_idx]);
                }
            }
            println!();
        }
    }
}