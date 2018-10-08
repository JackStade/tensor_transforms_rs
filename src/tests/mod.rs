use SymmetryObject;
use {Transform, TransformVector, Transformable};

#[test]
fn transform_math_add() {
    let (transform1, transform2) = get_test_transforms();
    let test1 = transform1.transform(&transform2);
    let test2 = transform1 + transform2;
    assert_eq!(&test1.transform[..], &test2.transform[..]);
}

#[test]
fn transform_math_invert() {
    let (transform1, transform2) = get_test_transforms();
    let (inv1, inv2) = (-transform1.clone(), -transform2.clone());
    let test1 = transform1.transform(&inv1);
    let test2 = transform2.transform(&inv2);
    let ident = Transform::identity(2);
    assert_eq!(&test1.transform[..], &ident.transform[..]);
    assert_eq!(&test2.transform[..], &ident.transform[..]);
}

#[test]
fn test_positive_group() {
    let s1 = SymmetryObject::new(2, vec![1, 1, 2, 2]).unwrap();
    let t1 = Transform::identity(2).get_transforms_positive(s1.get_transforms());
    for t in t1 {
        assert!(t.is_rigid_motion());
    }
    let s2 = SymmetryObject::new(3, vec![1, 1, 2, 2, 3, 3]).unwrap();
    let t2 = Transform::identity(3).get_transforms_positive(s2.get_transforms());
    for t in t2 {
        assert!(t.is_rigid_motion());
    }
    let s3 = SymmetryObject::new(3, vec![1, 1, 1, 1, 2, 2]).unwrap();
    let t3 = Transform::identity(3).get_transforms_positive(s3.get_transforms());
    for t in t3 {
        assert!(t.is_rigid_motion());
    }
    let s4 = SymmetryObject::new(3, vec![1, 2, 1, 2, 3, 3]).unwrap();
    let t4 = Transform::identity(3).get_transforms_positive(s4.get_transforms());
    for t in t4 {
        assert!(t.is_rigid_motion());
    }
}

fn get_test_transforms() -> (Transform, Transform) {
    let vec1 = TransformVector { dim: 0, sign: 1 };
    let vec2 = TransformVector { dim: 1, sign: -1 };
    // a flip y transform
    let transform1 = Transform {
        n: 2,
        transform: vec![vec1, vec2],
    };
    // a rotation
    let transform2 = Transform {
        n: 2,
        transform: vec![vec2, vec1],
    };
    (transform1, transform2)
}
