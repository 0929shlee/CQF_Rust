use std::io::{stdout, Write};
use rand::seq::SliceRandom;
use rand::thread_rng;
use crate::model::{Matrix, CqiMatrix, comp_quality};

pub struct ConnectionMatrix {
    pub cqi_matrix: CqiMatrix,
    pub matrix: Matrix,
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
        }
    }
    pub fn write(&self) {
        self.matrix.write();
    }

    pub fn brute_force(&mut self) {
        self.matrix.file_path = String::from("connection_matrix_brute_force.txt");
        for t in 0..self.matrix.n_time {
            println!("Brute Force Algorithm is running... {}%", (t + 1) * 100 / self.matrix.n_time);
            self._bf(t);
        }
        println!();
    }

    pub fn cqi_sorting(&mut self) {
        self.matrix.file_path = String::from("connection_matrix_cqi_sorting.txt");
        let mut coord_cqi_table = vec![vec![0 as u8; 4]; self.matrix.n_gnb * self.matrix.n_ue];

        for t in 0..self.matrix.n_time {
            println!("CQI Sorting Algorithm is running... {}%", (t + 1) * 100 / self.matrix.n_time);
            self._cs_get_coord_cqi_table(t, &mut coord_cqi_table);
            self._cs_connect_good_connection_candidates(t, &mut coord_cqi_table);
            self._cs_swap_connection(t);
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
            res += comp_quality::reward(&self.cqi_matrix, *g, u, time_idx);
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
            res += comp_quality::reward(&self.cqi_matrix, *g, u, time_idx);
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
    fn _cs_get_coord_cqi_table(&self, time_idx: usize, coord_cqi_table: &mut Vec<Vec<u8>>) {
        let mut cnt = 0;

        for g in 0..self.matrix.n_gnb {
            for u in 0..self.matrix.n_ue {
                coord_cqi_table[cnt][0] = g as u8;
                coord_cqi_table[cnt][1] = u as u8;
                coord_cqi_table[cnt][2] = self.cqi_matrix.matrix.matrix[g][u][time_idx];
                coord_cqi_table[cnt][3] = comp_quality::reward(
                    &self.cqi_matrix,
                    g,
                    u,
                    time_idx) as u8;
                cnt += 1;
            }
        }
    }
    fn _cs_connect_good_connection_candidates(
        &mut self,
        time_idx: usize,
        coord_cqi_table: &mut Vec<Vec<u8>>) {

        let mut gnb_connection_counts = vec![0 as usize; self.matrix.n_gnb];
        let mut ue_connection_counts = vec![0 as usize; self.matrix.n_ue];

        coord_cqi_table.sort_by(
            |v0, v1|
                v1.last().unwrap().cmp(v0.last().unwrap())
        );

        for v in coord_cqi_table {
            let g = v[0];
            let u = v[1];
            let cqi = v[2];

            //check if it is valid
            if cqi == 0 {
                break;
            }
            if gnb_connection_counts[g as usize] >= self.matrix.max_gnb_connection {
                continue;
            }
            if ue_connection_counts[u as usize] >= 1 {
                continue;
            }

            //valid
            self.matrix.matrix[g as usize][u as usize][time_idx] = 1;
            gnb_connection_counts[g as usize] += 1;
            ue_connection_counts[u as usize] += 1;
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
                            comp_quality::reward(&self.cqi_matrix, tmp_g, u, time_idx) +
                            comp_quality::reward(&self.cqi_matrix, g, tmp_u, time_idx) -
                            comp_quality::reward(&self.cqi_matrix, tmp_g, tmp_u, time_idx)
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