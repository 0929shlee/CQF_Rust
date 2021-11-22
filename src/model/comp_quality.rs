use crate::{ConnectionMatrix, CqiMatrix};

pub fn calculate(connection_matrix: &ConnectionMatrix) -> f64 {
    let mut res: f64 = 0.0;
    for t in 0..connection_matrix.matrix.n_time {
        let mut comp_quality_of_time: f64 = 0.0;
        for g in 0..connection_matrix.matrix.n_gnb {
            for u in 0..connection_matrix.matrix.n_ue {
                if connection_matrix.matrix.matrix[g][u][t] == 1 {
                    comp_quality_of_time += reward(
                        &connection_matrix.cqi_matrix,
                        g,
                        u,
                        t
                    );
                }
            }
        }
        comp_quality_of_time /= connection_matrix.matrix.n_ue as f64;
        res += comp_quality_of_time
    }
    res /= connection_matrix.matrix.n_time as f64;

    res
}
pub fn reward(
    cqi_matrix: &CqiMatrix,
    gnb_idx: usize,
    ue_idx: usize,
    time_idx: usize
) -> f64 {

    _get_reward_of_cqi(cqi_matrix.matrix.max_cqi) +
        _reward(cqi_matrix, gnb_idx, ue_idx, time_idx) -
        _penalty(cqi_matrix, gnb_idx, ue_idx, time_idx)
}

fn _reward(
    cqi_matrix: &CqiMatrix,
    gnb_idx: usize,
    ue_idx: usize,
    time_idx: usize
) -> f64 {

    _get_reward_of_cqi(cqi_matrix.matrix.matrix[gnb_idx][ue_idx][time_idx])
}
fn _penalty(
    cqi_matrix: &CqiMatrix,
    gnb_idx: usize,
    ue_idx: usize,
    time_idx: usize
) -> f64 {

    if time_idx > 0 {
        _get_penalty_of_cqi(
            cqi_matrix.matrix.max_cqi,
            cqi_matrix.matrix.matrix[gnb_idx][ue_idx][time_idx])
    }
    else {
        0.0
    }
}

fn _get_reward_of_cqi(cqi: u8) -> f64 {
    cqi as f64
}
fn _get_penalty_of_cqi(max_cqi: u8, cqi: u8) -> f64 {
    _get_reward_of_cqi(max_cqi) - _get_reward_of_cqi(cqi)
}