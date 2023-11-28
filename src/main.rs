use rand::{self, Rng};

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let mut rows: u32 = 24;
    let mut cols: u32 = 80;

    if let Some(rows_arg) = args.get(1) {
        if let Ok(rows_arg) = rows_arg.parse::<u32>() {
            rows = rows_arg;
        }
    }
    if let Some(cols_arg) = args.get(2) {
        if let Ok(cols_arg) = cols_arg.parse::<u32>() {
            cols = cols_arg;
        }
    }

    let smooth_version = 3;

    let mut smoothness: f64 = ((f64::from(rows) * f64::from(cols))
        / (match smooth_version {
            3 => 512.0,
            _ => 128.0,
        }))
        + 1.0;

    if let Some(smoothness_arg) = args.get(3) {
        if let Ok(smoothness_arg) = smoothness_arg.parse::<f64>() {
            smoothness = smoothness_arg;
        }
    }

    let chars: Vec<char> = " .`-',:_;~=^\"+/<>|()\\i%l{}!*I[]cr?stvx17afzenuCJLoj23Y45FT9PS6XbdkwGhVZpqAyEg&U#DKO$8@Hm0RBWNQM".chars().collect();

    let mut rng = rand::thread_rng();

    let original_map: Vec<Vec<f64>> = (0..rows)
        .map(|_| (0..cols).map(|_| rng.gen()).collect())
        .collect();

    let smoothed = match smooth_version {
        1 => smooth_v1(&original_map, smoothness as u32),
        2 => smooth_v2(&original_map, smoothness as u32),
        3 => smooth_v3(&original_map, smoothness),
        4 => smooth_v4(&original_map, smoothness),
        _ => clone_map(&original_map),
    };

    display(&smoothed, &chars);
}

fn smooth_v1(map: &Vec<Vec<f64>>, iterations: u32) -> Vec<Vec<f64>> {
    let mut old: Vec<Vec<f64>>;
    let mut new: Vec<Vec<f64>> = clone_map(&map);

    for _ in 0..iterations {
        old = clone_map(&new);

        for row in 0..old.len() {
            for col in 0..old[row].len() {
                new[row][col] = (old[wrap(row as i32 - 1, 0, old.len() as i32) as usize]
                    [wrap(col as i32 - 1, 0, old[row].len() as i32) as usize]
                    + old[wrap(row as i32 - 1, 0, old.len() as i32) as usize][col]
                    + old[wrap(row as i32 - 1, 0, old.len() as i32) as usize]
                        [wrap(col as i32 + 1, 0, old[row].len() as i32) as usize]
                    + old[row][wrap(col as i32 - 1, 0, old[row].len() as i32) as usize]
                    + old[row][col]
                    + old[row][wrap(col as i32 + 1, 0, old[row].len() as i32) as usize]
                    + old[wrap(row as i32 + 1, 0, old.len() as i32) as usize]
                        [wrap(col as i32 - 1, 0, old[row].len() as i32) as usize]
                    + old[wrap(row as i32 + 1, 0, old.len() as i32) as usize][col]
                    + old[wrap(row as i32 + 1, 0, old.len() as i32) as usize]
                        [wrap(col as i32 + 1, 0, old[row].len() as i32) as usize])
                    / 9.0;
            }
        }
    }

    new
}

fn smooth_v2(map: &Vec<Vec<f64>>, iterations: u32) -> Vec<Vec<f64>> {
    let corner = 1.0 / ((2.0_f64).sqrt() * 2.0);
    let edge = 1.0 / 2.0;
    let center = 1.0 / 1.0;
    let divisor = corner * 4.0 + edge * 4.0 + center;

    let mut old: Vec<Vec<f64>>;
    let mut new: Vec<Vec<f64>> = clone_map(&map);

    for _ in 0..iterations {
        old = clone_map(&new);

        for row in 0..old.len() {
            for col in 0..old[row].len() {
                new[row][col] = (old[wrap(row as i32 - 1, 0, old.len() as i32) as usize]
                    [wrap(col as i32 - 1, 0, old[row].len() as i32) as usize]
                    * corner
                    + old[wrap(row as i32 - 1, 0, old.len() as i32) as usize][col] * edge
                    + old[wrap(row as i32 - 1, 0, old.len() as i32) as usize]
                        [wrap(col as i32 + 1, 0, old[row].len() as i32) as usize]
                        * corner
                    + old[row][wrap(col as i32 - 1, 0, old[row].len() as i32) as usize] * edge
                    + old[row][col] * center
                    + old[row][wrap(col as i32 + 1, 0, old[row].len() as i32) as usize] * edge
                    + old[wrap(row as i32 + 1, 0, old.len() as i32) as usize]
                        [wrap(col as i32 - 1, 0, old[row].len() as i32) as usize]
                        * corner
                    + old[wrap(row as i32 + 1, 0, old.len() as i32) as usize][col] * edge
                    + old[wrap(row as i32 + 1, 0, old.len() as i32) as usize]
                        [wrap(col as i32 + 1, 0, old[row].len() as i32) as usize]
                        * corner)
                    / divisor;
            }
        }
    }

    new
}

fn smooth_v3(map: &Vec<Vec<f64>>, radius: f64) -> Vec<Vec<f64>> {
    let mut smoothed = clone_map(&map);

    let radius_ceil = radius.ceil() as i32;
    let radius_squared = (radius * radius).ceil() as i32;

    let weight_function = |center_x: f64, center_y: f64, point_x: f64, point_y: f64| -> f64 {
        let delta_x = (center_x - point_x).abs();
        let delta_y = (center_y - point_y).abs();
        let distance = (delta_x * delta_x + delta_y * delta_y).sqrt();
        let multiplier = radius - distance;
        multiplier.max(0.0)
    };

    let weight_map: Vec<Vec<f64>> = (0..radius.ceil() as u32)
        .map(|y| {
            (0..radius.ceil() as u32)
                .map(|x| weight_function(0.0, 0.0, x as f64, y as f64))
                .collect()
        })
        .collect();

    for (origin_y, origin_row) in map.iter().enumerate() {
        let y_start = origin_y as i32 - radius_ceil + 1;
        let y_end = origin_y as i32 + radius_ceil;
        for (origin_x, _) in origin_row.iter().enumerate() {
            let mut total = 0.0;
            let x_start = origin_x as i32 - radius_ceil + 1;
            let x_end = origin_x as i32 + radius_ceil;
            for y in y_start..y_end {
                let wrapped_y = wrap(y, 0, map.len() as i32) as usize;
                let y_diff = (origin_y as i32 - y).abs();
                for x in x_start..x_end {
                    let x_diff = (origin_x as i32 - x).abs();
                    if (x_diff * x_diff + y_diff * y_diff) <= radius_squared {
                        let wrapped_x = wrap(x, 0, origin_row.len() as i32) as usize;
                        let value = map[wrapped_y][wrapped_x];
                        let weight = weight_map[(origin_y as i32 - y).abs() as usize]
                            [(origin_x as i32 - x).abs() as usize];
                        total += value * weight;
                    }
                }
            }
            smoothed[origin_y][origin_x] = total;
        }
    }

    smoothed
}

fn smooth_v4(map: &Vec<Vec<f64>>, radius: f64) -> Vec<Vec<f64>> {
    let mut smoothed = clone_map(&map);

    let radius_ceil = radius.ceil() as i32;
    let radius_squared = (radius * radius).ceil() as i32;

    let weight_function = |center_x: f64, center_y: f64, point_x: f64, point_y: f64| -> f64 {
        let delta_x = (center_x - point_x).abs();
        let delta_y = (center_y - point_y).abs();
        let distance = (delta_x * delta_x + delta_y * delta_y).sqrt();
        let multiplier = radius - distance;
        multiplier.max(0.0)
    };

    let weight_map: Vec<Vec<f64>> = (0..radius.ceil() as u32)
        .map(|y| {
            (0..radius.ceil() as u32)
                .map(|x| weight_function(0.0, 0.0, x as f64, y as f64))
                .collect()
        })
        .collect();

    for (origin_y, origin_row) in map.iter().enumerate() {
        let y_start = origin_y as i32 - radius_ceil + 1;
        let y_end = origin_y as i32 + radius_ceil;
        for (origin_x, _) in origin_row.iter().enumerate() {
            let mut total = 0.0;
            let x_start = origin_x as i32 - radius_ceil + 1;
            let x_end = origin_x as i32 + radius_ceil;
            for y in y_start..y_end {
                let wrapped_y = wrap(y, 0, map.len() as i32) as usize;
                let y_diff = (origin_y as i32 - y).abs();
                for x in x_start..x_end {
                    let x_diff = (origin_x as i32 - x).abs();
                    if (x_diff * x_diff + y_diff * y_diff) <= radius_squared {
                        let wrapped_x = wrap(x, 0, origin_row.len() as i32) as usize;
                        let value = map[wrapped_y][wrapped_x];
                        let weight = weight_map[(origin_y as i32 - y).abs() as usize]
                            [(origin_x as i32 - x).abs() as usize];
                        total += value * weight;
                    }
                }
            }
            smoothed[origin_y][origin_x] = total;
        }
    }

    smoothed
}

fn display(display_map: &Vec<Vec<f64>>, chars: &Vec<char>) {
    let (min, max) = display_map
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |state, row| {
            let result = row
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |state, val| {
                    (val.min(state.0), val.max(state.1))
                });
            (result.0.min(state.0), result.1.max(state.1))
        });

    for row in display_map {
        for val in row {
            let char_index =
                map(*val, min, max, 0.0, f64::from((chars.len() - 1) as u32)).floor() as usize;
            print!("{}", chars[char_index]);
        }
        println!();
    }
}

fn map<
    T: std::ops::Add<Output = T>
        + std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + Copy,
>(
    x: T,
    in_min: T,
    in_max: T,
    out_min: T,
    out_max: T,
) -> T {
    (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
}

fn wrap<
    T: std::ops::Add<Output = T>
        + std::ops::Sub<Output = T>
        + std::ops::Rem<Output = T>
        + std::cmp::PartialOrd
        + Copy,
>(
    x: T,
    min: T,
    max: T,
) -> T {
    ((x % (max - min)) + (max - min)) % (max - min) + min
}

fn clone_map(map: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    map.iter().map(|x| x.clone()).collect()
}
