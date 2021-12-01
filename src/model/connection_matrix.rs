use std::io::{stdout, Write};
use rand::seq::SliceRandom;
use rand::thread_rng;
use crate::model::{Matrix, CqiMatrix, comp_quality};

pub struct ConnectionMatrix {
    pub cqi_matrix: CqiMatrix,
    pub matrix: Matrix,
    pub algorithm_name: String,
}

impl ConnectionMatrix {
    pub fn new(cqi_matrix: CqiMatrix) -> Self {
        let n_gnb = cqi_matrix.matrix.n_gnb;
        let n_ue = cqi_matrix.matrix.n_ue;
        let n_time = cqi_matrix.matrix.n_time;
        let max_gnb_connection = cqi_matrix.matrix.max_gnb_connection;

        ConnectionMatrix {
            cqi_matrix,
            matrix: Matrix::new(
                n_gnb,
                n_ue,
                n_time,
                String::from("connection_matrix.txt"),
                max_gnb_connection,
            ),
            algorithm_name: String::from(""),
        }
    }
    pub fn write(&self) {
        self.matrix.write();
    }

    fn set_algorithm_name(&mut self, name: &str) {
        self.algorithm_name = String::from(name);
    }
    fn set_connection_matrix_file_path(&mut self) {
        self.matrix.file_path = format!("connection_matrix_{}_{}.txt",
                                        self.algorithm_name,
                                        self.cqi_matrix.matrix.max_gnb_connection);
    }
    fn print_algorithm_start_msg(&self) {
        println!("{} Algorithm start", self.algorithm_name);
    }
    fn print_algorithm_running_msg(&self, time_idx: usize) {
        println!("{} Algorithm is running... ({}/{})", self.algorithm_name, time_idx + 1, self.matrix.n_time);
    }
    pub fn brute_force(&mut self) {
        self.set_algorithm_name("brute_force");
        self.set_connection_matrix_file_path();
        self.print_algorithm_start_msg();
        for t in 0..self.matrix.n_time {
            self.print_algorithm_running_msg(t);
            self._bf(t);
        }
        println!();
    }

    pub fn cqi_sorting(&mut self) {
        self.set_algorithm_name("cqi_sorting");
        self.set_connection_matrix_file_path();
        self.print_algorithm_start_msg();
        for t in 0..self.matrix.n_time {
            self.print_algorithm_running_msg(t);
            let mut coord_cqi_vector: Vec<(usize, usize, u8, f64)> = Vec::new();
            self._cs_get_coord_cqi_table(t, &mut coord_cqi_vector);
            self._cs_connect_good_connection_candidates(t, &mut coord_cqi_vector);
            self._cs_swap_connection(t);
        }
        println!();
    }

    pub fn expectation_sorting(&mut self) {
        self.set_algorithm_name("expectation_sorting");
        self.set_connection_matrix_file_path();
        self.print_algorithm_start_msg();
        for t in 0..self.matrix.n_time {
            self.print_algorithm_running_msg(t);
            let mut coord_cqi_vector: Vec<(usize, usize, u8, f64)> = Vec::new();
            self._es_get_coord_cqi_table(t, &mut coord_cqi_vector);
            self._es_connect_good_connection_candidates(t, &mut coord_cqi_vector);
            self._es_swap_connection(t);
        }
        println!();

    }

    pub fn naive_approach(&mut self) {
        self.set_algorithm_name("naive_approach");
        self.set_connection_matrix_file_path();
        self.print_algorithm_start_msg();
        for t in 0..self.matrix.n_time {
            self.print_algorithm_running_msg(t);
            self._na(t);
        }
        println!();
    }

    pub fn shuffle_split(&mut self, repeat: u64) {
        self.matrix.file_path = format!("connection_matrix_shuffle_split_{}.txt", repeat);
        for t in 0..self.matrix.n_time {
            println!("Shuffle-Split Algorithm is running... {}%", (t + 1) * 100 / self.matrix.n_time);
            self._ss(t, repeat);
        }
        println!();
    }
}

impl ConnectionMatrix {
    //Naive Approach
    fn _na(&mut self, time_idx: usize) {
        let mut gnb_connection_counts = vec![0 as usize; self.matrix.n_gnb];

        for u in 0..self.matrix.n_ue {
            let mut gnb_idx_cqi_arr: Vec<(usize, u8)> = vec![];
            for g in 0..self.matrix.n_gnb {
                gnb_idx_cqi_arr.push((g, self.cqi_matrix.matrix.matrix[g][u][time_idx]));
            }
            gnb_idx_cqi_arr.sort_by(
                |t0, t1|
                    t1.1.partial_cmp(&t0.1).unwrap()
            );
            for v in gnb_idx_cqi_arr {
                let g = v.0;
                if gnb_connection_counts[g] < self.matrix.max_gnb_connection {
                    assert_eq!(self.matrix.matrix[g][u][time_idx], 0);
                    self.matrix.matrix[g][u][time_idx] = 1;
                    gnb_connection_counts[g] += 1;
                    break;
                }
            }
        }
    }
}

impl ConnectionMatrix {
    //Expectation Sorting
    fn _es_get_coord_cqi_table(&self, time_idx: usize, coord_cqi_vector: &mut Vec<(usize, usize, u8, f64)>) {
        let mut sum_cqi_of_ue = vec![0 as u64; self.matrix.n_ue];
        for (u, s) in sum_cqi_of_ue.iter_mut().enumerate() {
            for g in 0..self.matrix.n_gnb {
                *s += (self.cqi_matrix.matrix.matrix[g][u][time_idx] as u64);
            }
        }
        let sum_cqi_of_ue = sum_cqi_of_ue;

        for g in 0..self.matrix.n_gnb {
            for u in 0..self.matrix.n_ue {
                coord_cqi_vector.push((g,
                                       u,
                                       self.cqi_matrix.matrix.matrix[g][u][time_idx],
                                       self._es_get_expectation(g, u, time_idx, &sum_cqi_of_ue)));
            }
        }
    }
    fn _es_get_expectation(
        &self,
        gnb_idx: usize,
        ue_idx: usize,
        time_idx: usize,
        sum_cqi_of_ue: &Vec<u64>
    ) -> f64 {
        comp_quality::_reward(
            &self.cqi_matrix,
            gnb_idx,
            ue_idx,
            time_idx
        )
            * self.matrix.n_gnb as f64
            - sum_cqi_of_ue[ue_idx] as f64
    }
    fn _es_connect_good_connection_candidates(
        &mut self,
        time_idx: usize,
        coord_cqi_vector: &mut Vec<(usize, usize, u8, f64)>) {

        let mut gnb_connection_counts = vec![0 as usize; self.matrix.n_gnb];
        let mut ue_connection_counts = vec![0 as usize; self.matrix.n_ue];

        coord_cqi_vector.sort_by(
            |t0, t1|
                t1.3.partial_cmp(&t0.3).unwrap()
        );

        for v in coord_cqi_vector {
            let g = v.0;
            let u = v.1;
            let cqi = v.2;

            //check if it is valid
            if cqi == 0 {
                break;
            }
            if gnb_connection_counts[g] >= self.matrix.max_gnb_connection {
                continue;
            }
            if ue_connection_counts[u] >= 1 {
                continue;
            }

            //valid
            self.matrix.matrix[g][u][time_idx] = 1;
            gnb_connection_counts[g] += 1;
            ue_connection_counts[u] += 1;
        }
    }
    fn _es_swap_connection(&mut self, time_idx: usize) {
        let mut gnb_connection_counts = vec![0 as usize; self.matrix.n_gnb];
        let mut ue_connection_counts = vec![0 as usize; self.matrix.n_ue];
        for g in 0..self.matrix.n_gnb {
            for u in 0..self.matrix.n_ue {
                if self.matrix.matrix[g][u][time_idx] == 1 {
                    gnb_connection_counts[g] += 1;
                    ue_connection_counts[u] += 1;
                }
            }
        }

        let mut sum_cqi_of_ue = vec![0 as u64; self.matrix.n_ue];
        for (u, s) in sum_cqi_of_ue.iter_mut().enumerate() {
            for g in 0..self.matrix.n_gnb {
                *s += (self.cqi_matrix.matrix.matrix[g][u][time_idx] as u64);
            }
        }

        for u in 0..self.matrix.n_ue {
            if ue_connection_counts[u] == 1 {
                continue;
            }
            let mut g: usize = 0;
            for c in gnb_connection_counts.iter() {
                if *c < self.matrix.max_gnb_connection {
                    break;
                }
                g += 1;
            }
            assert_eq!(self.matrix.matrix[g][u][time_idx], 0);

            if self.cqi_matrix.matrix.matrix[g][u][time_idx] != 0 {
                self.matrix.matrix[g][u][time_idx] = 1;
                continue;
            }

            let mut coord_cqi_vector: Vec<(usize, usize, f64)> = Vec::new(); //gnb_idx_to_go, ue_idx_to_switch, reward_changes
            for tmp_g in 0..self.matrix.n_gnb {
                if tmp_g == g {
                    continue;
                }
                for tmp_u in 0..self.matrix.n_ue {
                    if self.matrix.matrix[tmp_g][tmp_u][time_idx] == 1 {
                        coord_cqi_vector.push((
                            tmp_g,
                            tmp_u,
                            self._es_get_expectation(tmp_g, u, time_idx, &sum_cqi_of_ue) +
                                self._es_get_expectation(g, tmp_u, time_idx, &sum_cqi_of_ue) -
                                self._es_get_expectation(tmp_g, tmp_u, time_idx, &sum_cqi_of_ue)
                        ))
                    }
                }
            }
            coord_cqi_vector.sort_by(
                |t0, t1|
                    t1.2.partial_cmp(&t0.2).unwrap()
            );

            let res_g = coord_cqi_vector.first().unwrap().0;
            let res_u = coord_cqi_vector.first().unwrap().1;
            self.matrix.matrix[res_g][u][time_idx] = 1;
            self.matrix.matrix[res_g][res_u][time_idx] = 0;
            self.matrix.matrix[g][res_u][time_idx] = 1;

            gnb_connection_counts[g] += 1;
            ue_connection_counts[u] += 1;
        }
    }
}

impl ConnectionMatrix {
    //Shuffle-Split
    fn _ss(&mut self, time_idx: usize, repeat: u64) {
        let mut gnb_vector = vec![0 as usize; self.matrix.n_gnb * self.matrix.max_gnb_connection];
        self._ss_get_best_shuffled_gnb_vector(time_idx, repeat, &mut gnb_vector);
        self._ss_generate_connection(time_idx, &gnb_vector);
    }
    fn _ss_generate_connection(&mut self, time_idx: usize, gnb_vector: &Vec<usize>) {
        for (u, g) in gnb_vector[..self.matrix.n_ue].iter().enumerate() {
            self.matrix.matrix[*g][u][time_idx] = 1;
        }
    }
    fn _ss_is_valid_gnb_vector(&self, time_idx: usize, gnb_vector: &Vec<usize>) -> bool {
        for (u, g) in gnb_vector[..self.matrix.n_ue].iter().enumerate() {
            if self.cqi_matrix.matrix.matrix[*g][u][time_idx] == 0 {
                return false;
            }
        }
        true
    }
    fn _ss_shuffle_gnb_vector(&self, time_idx: usize, gnb_vector: &mut Vec<usize>) {
        loop {
            gnb_vector.shuffle(&mut thread_rng());
            if self._ss_is_valid_gnb_vector(time_idx, gnb_vector) {
                break;
            }
        }
    }
    fn _ss_get_best_shuffled_gnb_vector(
        &self,
        time_idx: usize,
        repeat: u64,
        gnb_vector: &mut Vec<usize>
    ) {
        self._ss_init_gnb_vector(gnb_vector);

        let mut max_gnb_vector = vec![0 as usize; self.matrix.n_ue];
        let mut max_comp_quality = -1.0;
        for _ in 0..repeat {
            self._ss_shuffle_gnb_vector(time_idx, gnb_vector);
            let comp_quality = self._ss_get_comp_quality_of_gnb_vector(gnb_vector, time_idx);
            if max_comp_quality < comp_quality {
                max_comp_quality = comp_quality;
                for u in 0..self.matrix.n_ue {
                    max_gnb_vector[u] = gnb_vector[u];
                }
            }
        }

        for (u, g) in max_gnb_vector.iter().enumerate() {
            gnb_vector[u] = *g;
        }
    }
    fn _ss_init_gnb_vector(&self, gnb_vector: &mut Vec<usize>) {
        for (i, n) in gnb_vector.iter_mut().enumerate() {
            *n = i / self.matrix.max_gnb_connection;
        }
    }
    fn _ss_get_comp_quality_of_gnb_vector(&self, gnb_vector: &Vec<usize>, time_idx: usize) -> f64 {
        let mut res = 0.0;
        for (u, g) in gnb_vector[..self.matrix.n_ue].iter().enumerate() {
            res += comp_quality::get_sum_of_comp_quality(&self.cqi_matrix, *g, u, time_idx);
        }

        res
    }
}

impl ConnectionMatrix {
    //Brute Force
    fn _bf(&mut self, time_idx: usize) {
        let mut gnb_vector = vec![0 as usize; self.matrix.n_ue];
        let mut max_gnb_vector = vec![0 as usize; self.matrix.n_ue];
        let mut max_comp_quality = -1.0;
        loop {
            if self._bf_is_valid_gnb_vector(&gnb_vector, time_idx) {
                let comp_quality = self._bf_get_comp_quality_of_gnb_vector(&gnb_vector, time_idx);
                if max_comp_quality < comp_quality {
                    max_comp_quality = comp_quality;
                    for u in 0..self.matrix.n_ue {
                        max_gnb_vector[u] = gnb_vector[u];
                    }
                }
            }

            if self._bf_is_last(&gnb_vector) {
                break;
            }
            self._bf_next_candidate(&mut gnb_vector);
        }

        for u in 0..self.matrix.n_ue {
            let g = max_gnb_vector[u];
            self.matrix.matrix[g][u][time_idx] = 1;
        }
    }
    fn _bf_next_candidate_r(&self, gnb_vector: &mut Vec<usize>, idx: usize) {
        if idx >= gnb_vector.iter().count() {
            return;
        }
        if gnb_vector[idx] == (self.matrix.n_gnb - 1) {
            gnb_vector[idx] = 0;
            self._bf_next_candidate_r(gnb_vector, idx - 1);
        }
        else {
            gnb_vector[idx] += 1;
        }
    }
    fn _bf_next_candidate(&self, gnb_vector: &mut Vec<usize>) {
        self._bf_next_candidate_r(gnb_vector, self.matrix.n_ue - 1);
    }
    fn _bf_get_comp_quality_of_gnb_vector(&self, gnb_vector: &Vec<usize>, time_idx: usize) -> f64 {
        let mut res = 0.0;
        for (u, g) in gnb_vector.iter().enumerate() {
            res += comp_quality::get_sum_of_comp_quality(&self.cqi_matrix, *g, u, time_idx);
        }

        res
    }
    fn _bf_is_last(&self, gnb_vector: &Vec<usize>) -> bool {
        let mut is_last = true;
        for v in gnb_vector {
            if *v != (self.matrix.n_gnb - 1) {
                is_last = false;
                break;
            }
        }

        is_last
    }
    fn _bf_is_valid_gnb_vector(&self, gnb_vector: &Vec<usize>, time_idx: usize) -> bool {
        let mut gnb_counts = vec![0 as usize; self.matrix.n_gnb];
        for u in 0..self.matrix.n_ue {
            let g = gnb_vector[u];
            gnb_counts[g] += 1;

            if self.cqi_matrix.matrix.matrix[g][u][time_idx] == 0 {
                return false;
            }
        }
        for n in gnb_counts {
            if n > self.matrix.max_gnb_connection {
                return false;
            }
        }

        true
    }
}
impl ConnectionMatrix {
    //CQI Sorting
    fn _cs_get_coord_cqi_table(&self, time_idx: usize, coord_cqi_vector: &mut Vec<(usize, usize, u8, f64)>) {
        for g in 0..self.matrix.n_gnb {
            for u in 0..self.matrix.n_ue {
                coord_cqi_vector.push((g,
                                      u,
                                      self.cqi_matrix.matrix.matrix[g][u][time_idx],
                                      comp_quality::get_sum_of_comp_quality(
                                          &self.cqi_matrix,
                                          g,
                                          u,
                                          time_idx)));
            }
        }
    }
    fn _cs_connect_good_connection_candidates(
        &mut self,
        time_idx: usize,
        coord_cqi_vector: &mut Vec<(usize, usize, u8, f64)>) {

        let mut gnb_connection_counts = vec![0 as usize; self.matrix.n_gnb];
        let mut ue_connection_counts = vec![0 as usize; self.matrix.n_ue];

        coord_cqi_vector.sort_by(
            |t0, t1|
                t1.3.partial_cmp(&t0.3).unwrap()
        );

        for v in coord_cqi_vector {
            let g = v.0;
            let u = v.1;
            let cqi = v.2;

            //check if it is valid
            if cqi == 0 {
                break;
            }
            if gnb_connection_counts[g] >= self.matrix.max_gnb_connection {
                continue;
            }
            if ue_connection_counts[u] >= 1 {
                continue;
            }

            //valid
            self.matrix.matrix[g][u][time_idx] = 1;
            gnb_connection_counts[g] += 1;
            ue_connection_counts[u] += 1;
        }
    }
    fn _cs_swap_connection(&mut self, time_idx: usize) {
        let mut gnb_connection_counts = vec![0 as usize; self.matrix.n_gnb];
        let mut ue_connection_counts = vec![0 as usize; self.matrix.n_ue];
        for g in 0..self.matrix.n_gnb {
            for u in 0..self.matrix.n_ue {
                if self.matrix.matrix[g][u][time_idx] == 1 {
                    gnb_connection_counts[g] += 1;
                    ue_connection_counts[u] += 1;
                }
            }
        }

        for u in 0..self.matrix.n_ue {
            if ue_connection_counts[u] == 1 {
                continue;
            }
            let mut g: usize = 0;
            for c in gnb_connection_counts.iter() {
                if *c < self.matrix.max_gnb_connection {
                    break;
                }
                g += 1;
            }
            assert_eq!(self.matrix.matrix[g][u][time_idx], 0);

            if self.cqi_matrix.matrix.matrix[g][u][time_idx] != 0 {
                self.matrix.matrix[g][u][time_idx] = 1;
                continue;
            }

            let mut coord_cqi_matrix: Vec<(usize, usize, f64)> = Vec::new(); //gnb_idx_to_go, ue_idx_to_switch, reward_changes
            for tmp_g in 0..self.matrix.n_gnb {
                if tmp_g == g {
                    continue;
                }
                for tmp_u in 0..self.matrix.n_ue {
                    if self.matrix.matrix[tmp_g][tmp_u][time_idx] == 1 {
                        coord_cqi_matrix.push((
                            tmp_g,
                            tmp_u,
                            comp_quality::get_sum_of_comp_quality(&self.cqi_matrix, tmp_g, u, time_idx) +
                            comp_quality::get_sum_of_comp_quality(&self.cqi_matrix, g, tmp_u, time_idx) -
                            comp_quality::get_sum_of_comp_quality(&self.cqi_matrix, tmp_g, tmp_u, time_idx)
                        ))
                    }
                }
            }
            coord_cqi_matrix.sort_by(
                |t0, t1|
                    t1.2.partial_cmp(&t0.2).unwrap()
            );

            let res_g = coord_cqi_matrix.first().unwrap().0;
            let res_u = coord_cqi_matrix.first().unwrap().1;
            self.matrix.matrix[res_g][u][time_idx] = 1;
            self.matrix.matrix[res_g][res_u][time_idx] = 0;
            self.matrix.matrix[g][res_u][time_idx] = 1;

            gnb_connection_counts[g] += 1;
            ue_connection_counts[u] += 1;
        }
    }
}

impl ConnectionMatrix {
    //test
    pub fn is_valid(&self) -> bool {
        let mut is_valid = true;

        for t in 0..self.matrix.n_time {
            let mut sum_of_ue_cqi_vector = vec![0 as u8; self.matrix.n_ue];
            let mut n_ue_connected_to_gnb_vector = vec![0 as u8; self.matrix.n_gnb];
            let mut n_gnb_connected_to_ue_vector = vec![0 as u8; self.matrix.n_ue];

            for u in 0..self.matrix.n_ue {
                for g in 0..self.matrix.n_gnb {
                    sum_of_ue_cqi_vector[u] += self.cqi_matrix.matrix.matrix[g][u][t];
                    n_gnb_connected_to_ue_vector[u] += self.matrix.matrix[g][u][t];
                    n_ue_connected_to_gnb_vector[g] += self.matrix.matrix[g][u][t];
                }
            }

            let n_ue_zero_cqi = sum_of_ue_cqi_vector
                .iter()
                .filter(|n| **n == 0)
                .count();
            let n_ue_connected = n_gnb_connected_to_ue_vector
                .iter()
                .filter(|n| **n >= 1)
                .count();

            if n_ue_connected_to_gnb_vector
                .iter()
                .filter(|n| **n > self.matrix.max_gnb_connection as u8)
                .count() > 0
                ||
                n_ue_zero_cqi + n_ue_connected < self.matrix.n_ue {

                is_valid = false;
                println!("\nerror: Invalid connection detected at time {}.", t + 1);
                println!("CQI Matrix: ");
                self.cqi_matrix.matrix.print_of_time(t);
                println!("Connection Matrix: ");
                self.matrix.print_of_time(t);
                println!();
            }
        }

        is_valid
    }
}