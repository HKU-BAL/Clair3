use std::borrow::Borrow;

pub(crate) trait KhKey: Eq {
    fn hash(&self) -> u32;
}

impl KhKey for String {
    fn hash(&self) -> u32 {
        self.as_bytes().iter().fold(0_u32, |hash, byte| {
            hash.wrapping_mul(31).wrapping_add(u32::from(*byte))
        })
    }
}

impl KhKey for usize {
    fn hash(&self) -> u32 {
        *self as u32
    }
}

impl KhKey for str {
    fn hash(&self) -> u32 {
        self.as_bytes().iter().fold(0_u32, |hash, byte| {
            hash.wrapping_mul(31).wrapping_add(u32::from(*byte))
        })
    }
}

#[derive(Debug, Clone)]
pub(crate) struct KhCounter<K> {
    buckets: Vec<Option<(K, i64)>>,
    occupied: usize,
    upper_bound: usize,
}

impl<K> Default for KhCounter<K> {
    fn default() -> Self {
        Self {
            buckets: Vec::new(),
            occupied: 0,
            upper_bound: 0,
        }
    }
}

impl<K: KhKey> KhCounter<K> {
    fn upper_bound(size: usize) -> usize {
        (size as f64 * 0.77 + 0.5) as usize
    }

    fn resize(&mut self, requested: usize) {
        let new_size = requested.max(4).next_power_of_two();
        let mut old = std::mem::take(&mut self.buckets);
        let old_size = old.len();
        let new_mask = new_size - 1;
        let mut moved = Vec::with_capacity(new_size);
        moved.resize_with(new_size, || None);
        let mut occupied = vec![false; new_size];
        for old_index in 0..old_size {
            let Some(mut entry) = old[old_index].take() else {
                continue;
            };
            loop {
                let mut slot = entry.0.hash() as usize & new_mask;
                let mut step = 0;
                while occupied[slot] {
                    step += 1;
                    slot = (slot + step) & new_mask;
                }
                occupied[slot] = true;
                if slot < old_size && old[slot].is_some() {
                    let mut previous = old[slot].take().unwrap();
                    std::mem::swap(&mut previous, &mut entry);
                    moved[slot] = Some(previous);
                } else {
                    moved[slot] = Some(entry);
                    break;
                }
            }
        }
        self.occupied = moved.iter().filter(|entry| entry.is_some()).count();
        self.buckets = moved;
        self.upper_bound = Self::upper_bound(new_size);
    }

    pub fn entry(&mut self, key: K) -> &mut i64 {
        if self.occupied >= self.upper_bound {
            let requested = if self.buckets.len() > self.occupied.saturating_mul(2) {
                self.buckets.len() - 1
            } else {
                self.buckets.len() + 1
            };
            self.resize(requested);
        }
        let mask = self.buckets.len() - 1;
        let mut index = key.hash() as usize & mask;
        let mut step = 0;
        loop {
            match &self.buckets[index] {
                Some((existing, _)) if existing == &key => {
                    return &mut self.buckets[index].as_mut().unwrap().1;
                }
                Some(_) => {
                    step += 1;
                    index = (index + step) & mask;
                }
                None => {
                    self.buckets[index] = Some((key, 0));
                    self.occupied += 1;
                    return &mut self.buckets[index].as_mut().unwrap().1;
                }
            }
        }
    }

    pub fn get<Q>(&self, key: &Q) -> Option<&i64>
    where
        K: Borrow<Q>,
        Q: KhKey + ?Sized,
    {
        if self.buckets.is_empty() {
            return None;
        }
        let mask = self.buckets.len() - 1;
        let mut index = key.hash() as usize & mask;
        let mut step = 0;
        loop {
            match &self.buckets[index] {
                Some((existing, value)) if existing.borrow() == key => return Some(value),
                Some(_) => {
                    step += 1;
                    index = (index + step) & mask;
                }
                None => return None,
            }
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = (&K, &i64)> {
        self.buckets
            .iter()
            .filter_map(|entry| entry.as_ref().map(|(key, value)| (key, value)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hashes_and_resize_preserves_counter_values() {
        let mut counter = KhCounter::default();
        for key in ["TA", "TAAAAA", "TAA", "G", "CT", "GA", "CAA", "CAAAAAA"] {
            *counter.entry(key.to_string()) += 1;
        }
        for key in ["TA", "TAAAAA", "TAA", "G", "CT", "GA", "CAA", "CAAAAAA"] {
            assert_eq!(counter.get(&key.to_string()), Some(&1));
        }
        assert_eq!(counter.iter().count(), 8);
        assert_eq!(
            counter
                .iter()
                .map(|(key, _)| key.clone())
                .collect::<Vec<_>>(),
            ["CT", "CAA", "TAA", "CAAAAAA", "G", "GA", "TA", "TAAAAA"]
        );
        let mut integers = KhCounter::default();
        for key in 1..=12 {
            *integers.entry(key) += 1;
        }
        assert_eq!(integers.get(&12), Some(&1));
        assert_eq!(
            integers.iter().map(|(key, _)| *key).collect::<Vec<_>>(),
            (1..=12).collect::<Vec<_>>()
        );
    }

    #[test]
    fn repeated_insertions_survive_multiple_resizes() {
        let mut counter = KhCounter::default();
        for (key, count) in [
            ("TAA", 10),
            ("TAAAAA", 13),
            ("TA", 3),
            ("TAAA", 3),
            ("TAAAAAAAA", 1),
            ("TAAAAAA", 4),
            ("TAAAAAAAAA", 2),
        ] {
            for _ in 0..count {
                *counter.entry(key.to_string()) += 1;
            }
            assert_eq!(
                counter.get(&key.to_string()),
                Some(&count),
                "count for {key}"
            );
        }
        assert_eq!(
            counter
                .iter()
                .map(|(key, _)| key.clone())
                .collect::<Vec<_>>(),
            [
                "TA",
                "TAAAAAAAAA",
                "TAAAAAA",
                "TAA",
                "TAAAAAAAA",
                "TAAA",
                "TAAAAA"
            ]
        );
    }
}
