use crate::model::{Matrix};
use rand::Rng;

pub struct CqiMatrix {
    scaling_factor: f64,
    noise_density_factor: f64,
    noise_length_factor: usize,
    pub matrix: Matrix,
}

impl CqiMatrix {
    pub fn new(
        n_gnb: usize,
        n_ue: usize,
        n_time: usize,
        max_gnb_connection: usize,
        scaling_factor: f64,
        noise_density_factor: f64,
        noise_length_factor: usize) -> Self {

        CqiMatrix {
            scaling_factor,
            noise_density_factor,
            noise_length_factor,
            matrix: Matrix::new(
                n_gnb,
                n_ue,
                n_time,
                String::from("cqi_matrix.txt"),
                max_gnb_connection,
            ),
        }
    }
    pub fn write(&self) {
        self.matrix.write();
    }

    pub fn generate(&mut self) {
        if self.matrix.is_there_file() {
            self.matrix.read();
        }
        else {
            self.generate_random_matrix();
            self.write();
        }
    }
    fn generate_random_matrix(&mut self) {
        for g in 0..self.matrix.n_gnb {
            for u in 0..self.matrix.n_ue {
                let n_turning_pos = self.generate_n_turning_pos();

                let mut turning_values: Vec<u8> = vec![0; n_turning_pos + 2];
                self.generate_n_turning_values(&mut turning_values);

                let mut turning_pos: Vec<usize> = vec![0; n_turning_pos + 2];
                self.generate_turning_pos(&mut turning_pos);

                let mut cqi_vector: Vec<u8> = vec![0; self.matrix.n_time];
                self.generate_cqi_vector(
                    &turning_values,
                    &turning_pos,
                    &mut cqi_vector);

                self.add_noise_to_cqi_vector(&mut cqi_vector);

                for t in 0..self.matrix.n_time {
                    self.matrix.matrix[g][u][t] = cqi_vector[t];
                }
            }
        }
    }
}

impl CqiMatrix {
    fn generate_n_turning_pos(&self) -> usize {
        let mut rng = rand::thread_rng();
        let range = self.matrix.n_time as f64 * self.scaling_factor - 1.0;
        let range = range as usize;
        let res: usize = rng.gen_range(0..range);

        res
    }
    fn generate_n_turning_values(&self, turning_values: &mut Vec<u8>) {
        let mut rng = rand::thread_rng();
        for n in turning_values {
            *n = rng.gen_range(0..(self.matrix.max_cqi + 1)) as u8;
        }
    }
    fn generate_turning_pos(&self, turning_pos: &mut Vec<usize>) {
        let length = turning_pos.iter().count();
        for i in 0..(length - 1) {
            turning_pos[i] = (self.matrix.n_time / (length - 1)) * i;
        }
        turning_pos[length - 1] = self.matrix.n_time - 1;
    }
    fn generate_cqi_vector(
        &self,
        turning_values: &Vec<u8>,
        turning_pos: &Vec<usize>,
        cqi_vector: &mut Vec<u8>) {

        let length = turning_pos.iter().count();
        for i in 0..(length - 1) {
            let dist = turning_pos[i + 1] - turning_pos[i];
            if turning_values[i] <= turning_values[i + 1] { // increasing
                let diff = (turning_values[i + 1] - turning_values[i]) as usize;
                for j in 0..dist {
                    cqi_vector[turning_pos[i] + j] = turning_values[i] + (diff * j / dist) as u8;
                }
            }
            else { // decreasing
                let diff = (turning_values[i] - turning_values[i + 1]) as usize;
                for j in 0..dist {
                    cqi_vector[turning_pos[i] + j] = turning_values[i] - (diff * j / dist) as u8;
                }
            }
        }
        *(cqi_vector.last_mut().unwrap()) = *(turning_values.last().unwrap());
    }
    fn add_noise_to_cqi_vector(&self, cqi_vector: &mut Vec<u8>) {
        let range = (self.matrix.n_time as f64 * self.noise_density_factor) as usize;
        let mut rng = rand::thread_rng();
        for _ in 0..range {
            let rand_pos = rng.gen_range(0..self.matrix.n_time);
            for j in 0..self.noise_length_factor {
                let rand_value = rng.gen_range(0..(cqi_vector[rand_pos + j] + 1));
                cqi_vector[rand_pos + j] -= rand_value;

                if &(cqi_vector[rand_pos + j]) == cqi_vector.last().unwrap() {
                    break;
                }
            }
        }
    }

}