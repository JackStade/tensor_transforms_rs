use {Transform, TransformDimension, Transformable};

/// Simmilar to a `SymmetryObject`, this type has one value for every direction.
///
/// A DirectionObject cannot be transformed by a transform with a higher number of dimensions,
/// because it would be unable to fill empty transforms with a default value. In the future, there might
/// be a seperate implementation for `DirectionObject`s where T: Default, however for the moment this
/// is not easily possible with Rust.
#[derive(Clone)]
pub struct DirectionObject<T> {
    n: usize,
    vals: Vec<T>,
}

impl<U> DirectionObject<U> {
    pub fn new<T>(vals: Vec<T>) -> Result<DirectionObject<T>, &'static str> {
        if (&vals).len() % 2 != 0 {
            return Err("Length must be a multiple of 2.");
        }
        Ok(DirectionObject {
            n: (&vals).len() / 2,
            vals: vals,
        })
    }

    pub fn get_vals(&self) -> &[U] {
        &self.vals[..]
    }
}

impl<T> TransformDimension for DirectionObject<T> {
    #[inline(always)]
    fn dimensions(&self) -> usize {
        self.n
    }
}

impl<T: PartialEq + Sized + Clone> Transformable for DirectionObject<T> {
    fn transform(&self, transform: &Transform) -> DirectionObject<T> {
        let n = transform.dimensions();
        if n > self.n {
            panic!("Cannot transform a direction object into a higher dimension")
        };
        let axes = transform.values();
        let mut new_vals = Vec::with_capacity(2 * self.n);
        unsafe {
            new_vals.set_len(2 * self.n);
        }
        for i in 0..self.n {
            let (dim, dir) = if i < n {
                (2 * axes[i].dim, ((axes[i].sign + 1) / 2) as usize)
            } else {
                (2 * i, 0)
            };
            new_vals[dim + dir] = self.vals[dim].clone();
            new_vals[dim + 1 - dir] = self.vals[dim + 1].clone();
        }
        DirectionObject {
            n: self.n,
            vals: new_vals,
        }
    }
}

/// A simple n-dimensional tensor that can be transformed.
///
/// This type can be transformed by a transform of any size,
/// since it does not need to add extra values when the dimension is increased.
///
/// For the moment, this type only implements transform if T is `Clone`. In the future
/// there might be a way to transform it in place.
#[derive(Clone)]
pub struct TransformTensor<T: Sized> {
    // number of dimensions
    n: usize,
    // the size of each dimension
    size: Vec<usize>,
    // values
    vals: Vec<T>,
}

impl<T: Sized> TransformTensor<T> {
    /// Generates a new `TransformTensor`, checking to make sure the Vec has enough elements.
    pub fn new(size: &[usize], vals: Vec<T>) -> Result<TransformTensor<T>, &'static str> {
        let mut product = 1;
        for i in size {
            product *= i
        }
        if Vec::len(&vals) != product {
            return Err("Vector does not have enough elements.");
        }
        Ok(TransformTensor {
            n: size.len(),
            size: size.to_vec(),
            vals: vals,
        })
    }

    pub fn get_size(&self) -> &[usize] {
        &self.size[..]
    }

    pub fn get_vals(&self) -> &[T] {
        &self.vals[..]
    }
}

use std::cmp::max;

impl<T> TransformDimension for TransformTensor<T> {
    #[inline(always)]
    fn dimensions(&self) -> usize {
        self.n
    }
}

impl<T: Sized + Clone> Transformable for TransformTensor<T> {
    fn transform(&self, transform: &Transform) -> Self {
        let n = transform.dimensions();
        let max = max(n, self.n);
        let axes = transform.values();
        let mut new_size = vec![1; max];
        for i in 0..self.n {
            let sindex = if i < n { axes[i].dim } else { i };
            new_size[sindex] = self.size[i];
        }
        let mut new_page_sizes = Vec::with_capacity(max + 1);
        let mut product = 1;
        new_page_sizes.push(product);
        for i in 0..max {
            product *= new_size[i];
            new_page_sizes.push(product);
        }
        let mut product = 1;
        let mut old_page_sizes = Vec::with_capacity(self.n + 1);
        old_page_sizes.push(product);
        for i in 0..max {
            product *= self.size[i];
            old_page_sizes.push(product);
        }

        let len = Vec::len(&self.vals);
        let mut new_vals = Vec::with_capacity(len);
        // used to keep track of which elements are initialzed during debugging
        let mut num_set = 0;
        let mut vals_set = if cfg!(debug_assertions) {
            Vec::with_capacity(len)
        } else {
            Vec::new()
        };
        // set the length of new_vals
        unsafe {
            new_vals.set_len(len);
        }
        let mut loop_indices = vec![0; self.n];
        let mut offset = vec![0; self.n];
        let mut transform_coordinate = 0;
        for i in 0..self.n {
            let (dim, dir) = if i < n {
                (axes[i].dim, axes[i].sign)
            } else {
                (i, 1)
            };
            offset[i] = new_page_sizes[dim] as isize * dir as isize;
            // start at the opposite end for negative indices
            if dir == -1 {
                transform_coordinate += (new_page_sizes[dim] * (self.size[i] - 1)) as isize;
            }
        }

        let mut i = 0;
        let mut old_index = 0;
        while i < self.n {
            if loop_indices[i] < self.size[i] {
                if i == 0 {
                    // if there are debug assertions, we check to make sure that each element of the
                    // vector is set once and only once.
                    if cfg!(debug_assertions) {
                        if vals_set[transform_coordinate as usize] {
                            panic!(
                                "DEBUG: Transform set value: {} more than once.",
                                transform_coordinate
                            );
                        } else {
                            vals_set[transform_coordinate as usize] = true;
                            num_set += 1;
                        }
                    }
                    new_vals[transform_coordinate as usize] = self.vals[old_index].clone();
                    old_index += 1;
                }
                loop_indices[i] += 1;

                transform_coordinate += offset[i];
                if loop_indices[i] < self.size[i] {
                    i = 0
                };
            } else {
                transform_coordinate -= offset[i] * self.size[i] as isize;
                loop_indices[i] = 0;
                i += 1;
            }
        }
        // check that the total number of elements set is num
        if cfg!(debug_assertions) && num_set != len {
            panic!(
                "DEBUG: Transform set only {} values out of {}.",
                num_set, len
            );
        }
        TransformTensor {
            n: max,
            size: new_size,
            vals: new_vals,
        }
    }
}

use fmt::{Display, Formatter};
use std::fmt;

impl<T: Sized + Display> Display for TransformTensor<T> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        if self.n == 1 {
            writeline(&self.vals[..], f)?;
        } else if self.n == 2 {
            writelines(&self.vals[..], self.size[0], f)?;
        } else {
            let page_size = self.size[0] * self.size[1];
            let line_size = self.size[0];
            let mut product = 1;
            for i in 2..self.n {
                product *= self.size[i];
            }
            for i in 0..product {
                writelines(&self.vals[page_size * i..page_size * (i + 1)], line_size, f)?;
                if i < product - 1 {
                    write!(f, "\n\n")?
                };
            }
        }
        Ok(())
    }
}

fn writelines<T: Display>(vals: &[T], line_length: usize, f: &mut Formatter) -> fmt::Result {
    let num_lines = vals.len() / line_length;
    for i in 0..num_lines {
        writeline(&vals[i * line_length..(i + 1) * line_length], f)?;
        if i < num_lines - 1 {
            write!(f, "\n")?
        };
    }
    Ok(())
}

fn writeline<T: Display>(vals: &[T], f: &mut Formatter) -> fmt::Result {
    for val in &vals[..vals.len() - 1] {
        write!(f, "{}, ", val)?;
    }
    write!(f, "{}", vals[vals.len() - 1])?;
    Ok(())
}
