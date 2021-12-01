use crate::{ConnectionMatrix, CqiMatrix};

pub fn calculate(connection_matrix: &ConnectionMatrix) -> f64 {
    let mut res: f64 = 0.0;
    for t in 0..connection_matrix.matrix.n_time {
        let mut comp_quality_of_time: f64 = 0.0;
        for g in 0..connection_matrix.matrix.n_gnb {
            for u in 0..connection_matrix.matrix.n_ue {
                if connection_matrix.matrix.matrix[g][u][t] == 1 {
                    comp_quality_of_time += get_sum_of_comp_quality(
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
pub fn get_sum_of_comp_quality(
    cqi_matrix: &CqiMatrix,
    gnb_idx: usize,
    ue_idx: usize,
    time_idx: usize
) -> f64 {
    _max_reward(cqi_matrix) +
        _reward(cqi_matrix, gnb_idx, ue_idx, time_idx) -
        _penalty(cqi_matrix, gnb_idx, ue_idx, time_idx)
}

pub fn _max_reward(
    cqi_matrix: &CqiMatrix,
) -> f64 {

    _get_reward_of_cqi(cqi_matrix.matrix.max_cqi)
}
pub fn _reward(
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
        cqi_matrix.penalty_factor * (_get_reward_of_cqi(cqi_matrix.matrix.max_cqi) -
            _get_reward_of_cqi(cqi_matrix.matrix.matrix[gnb_idx][ue_idx][time_idx - 1]))
    }
    else {
        0.0
    }
}

fn _get_reward_of_cqi(cqi: u8) -> f64 {
    /*
    let arr = [
        0.0,
        0.2137962089502232,
        0.33884415613920255,
        0.588843655355589,
        1.0471285480508996,
        1.7378008287493754,
        2.6915348039269156,
        3.8904514499428067,
        6.456542290346555,
        10.715193052376065,
        14.791083881682072,
        25.703957827688633,
        42.657951880159274,
        74.13102413009173,
        125.89254117941675,
        186.20871366628674
    ];
    arr[cqi as usize]

     */
    /*
    same as
     */
    get_sinr_factor(cqi_to_sinr(cqi))
}

fn get_sinr_factor(sinr: f64) -> f64 {
    10.0_f64.powf(sinr / 10.0)
}
fn cqi_to_sinr(cqi: u8) -> f64 {
    let sinr_arr = [-100.0, -6.7, -4.7, -2.3, 0.2, 2.4, 4.3, 5.9, 8.1, 10.3, 11.7, 14.1, 16.3, 18.7, 21.0, 22.7];
    assert!((cqi as usize) < sinr_arr.len());
    sinr_arr[cqi as usize]
}
