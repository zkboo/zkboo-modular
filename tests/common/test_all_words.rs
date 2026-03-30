#[macro_export]
#[doc(hidden)]
macro_rules! _on_given_words {
    (
        [$($W: ty),* $(,)?],
        $m:ident ! ( $($args:tt)* )
    ) => {
        ::paste::paste! {
            $m!(
                {$([<$W>] : ($W, 1, $W)),*,},
                $($args)*
            );
        }
    };
}

#[doc(inline)]
pub use _on_given_words as on_given_words;

#[macro_export]
#[doc(hidden)]
macro_rules! _on_all_words {
    (
        $m:ident ! ( $($args:tt)* )
    ) => {
        ::paste::paste! {
            $crate::common::test_all_words::on_given_words!(
                [u8, u16, u32, u64, u128],
                $m ! ($($args)* )
            );
        }
    };
}

#[doc(inline)]
pub use _on_all_words as on_all_words;

#[macro_export]
#[doc(hidden)]
macro_rules! _test_on_all_words {
    (
        $func: ident,
        $m:ident ! ( $($args:tt)* )
    ) => {
    ::paste::paste! {
            #[test]
            fn [<test_ $func>]() {
                $crate::common::test_all_words::on_all_words!(
                    $m ! ([<test_ $func>],  $($args)* )
                );
            }
        }
    };
}

#[doc(inline)]
pub use _test_on_all_words as test_on_all_words;

#[macro_export]
#[doc(hidden)]
macro_rules! _on_given_composites {
    (
        [$($W: ty),* $(,)?],
        $m:ident ! ( $($args:tt)* )
    ) => {
        ::paste::paste! {
            $m!(
                {
                    $([<$W _ 1>] : ($W, 1, ::zkboo::word::CompositeWord::<$W, 1>)),*,
                    $([<$W _ 2>] : ($W, 2, ::zkboo::word::CompositeWord::<$W, 2>)),*,
                    $([<$W _ 3>] : ($W, 3, ::zkboo::word::CompositeWord::<$W, 3>)),*,
                    $([<$W _ 4>] : ($W, 4, ::zkboo::word::CompositeWord::<$W, 4>)),*,
                },
                $($args)*
            );
        }
    };
}

#[doc(inline)]
pub use _on_given_composites as on_given_composites;

#[macro_export]
#[doc(hidden)]
macro_rules! _on_all_composites {
    (
        $m:ident ! ( $($args:tt)* )
    ) => {
        ::paste::paste! {
            $crate::common::test_all_words::on_given_composites!(
                [u8, u16, u32, u64, u128],
                $m ! ($($args)* )
            );
        }
    };
}

#[doc(inline)]
pub use _on_all_composites as on_all_composites;

#[macro_export]
#[doc(hidden)]
macro_rules! _on_all_words_and_composites {
    (
        $m:ident ! ( $($args:tt)* )
    ) => {
        ::paste::paste! {
            $crate::common::test_all_words::on_given_words!(
                [u8, u16, u32, u64, u128],
                $m ! ($($args)* )
            );
            $crate::common::test_all_words::on_given_composites!(
                [u8, u16, u32, u64, u128],
                $m ! ($($args)* )
            );
        }
    };
}

#[doc(inline)]
pub use _on_all_words_and_composites as on_all_words_and_composites;

#[macro_export]
#[doc(hidden)]
macro_rules! _test_on_all_composites {
    (
        $func: ident,
        $m:ident ! ( $($args:tt)* )
    ) => {
    ::paste::paste! {
            #[test]
            fn [<test_ $func>]() {
                $crate::common::test_all_words::on_all_composites!(
                    $m ! ([<test_ $func>],  $($args)* )
                );
            }
        }
    };
}

#[doc(inline)]
pub use _test_on_all_composites as test_on_all_composites;

#[macro_export]
#[doc(hidden)]
macro_rules! _test_on_all_words_and_composites {
    (
        $func: ident,
        $m:ident ! ( $($args:tt)* )
    ) => {
    ::paste::paste! {
            #[test]
            fn [<test_ $func>]() {
                $crate::common::test_all_words::on_all_words!(
                    $m ! ([<test_ $func>],  $($args)* )
                );
                $crate::common::test_all_words::on_all_composites!(
                    $m ! ([<test_ $func>],  $($args)* )
                );
            }
        }
    };
}

#[doc(inline)]
pub use _test_on_all_words_and_composites as test_on_all_words_and_composites;
