use rand::thread_rng;
use rand::Rng;

// Partition array around pivot, return pivot's final index
pub fn partition(arr: &mut [f64], pivot_index: usize) -> usize {
  arr.swap(pivot_index, arr.len() - 1);
  let mut store_index = 0;
  for i in 0..arr.len() - 1 {
    if arr[i] < arr[arr.len() - 1] {
      arr.swap(store_index, i);
      store_index += 1;
    }
  }
  arr.swap(store_index, arr.len() - 1);
  store_index
}

// Find k-th smallest element using Quickselect
pub fn quickselect(arr: &mut [f64], k: usize) -> f64 {
  if arr.len() == 1 {
    return arr[0];
  }
  let pivot_index = thread_rng().gen_range(0..arr.len());
  let pivot_index = partition(arr, pivot_index);

  match pivot_index.cmp(&k) {
    std::cmp::Ordering::Equal => arr[pivot_index],
    std::cmp::Ordering::Less => quickselect(&mut arr[pivot_index + 1..], k - pivot_index - 1),
    std::cmp::Ordering::Greater => quickselect(&mut arr[..pivot_index], k),
  }
}

// Compute a percentile value from the data
pub fn percentile(data: &mut Vec<f64>, percentile: f64) -> f64 {
  assert!(percentile >= 0.0 && percentile <= 100.0);
  let k = (percentile / 100.0 * (data.len() - 1) as f64).round() as usize;
  quickselect(data, k)
}

