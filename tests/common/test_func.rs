#[macro_export]
#[doc(hidden)]
macro_rules! _define_test_circuit {
    (
        $name: ident,
        {$($in:ident : $in_t:ty),* $(,)?},
        [$($out:ident),* $(,)?],
        |$executor: ident| $body: block
        $(, $param_name: ident : $param_type: ty)?
    ) => {
        ::paste::paste! {
            struct $name{
                $(
                    $in: $in_t,
                )*
                $($param_name: $param_type),*
            }
            impl Default for $name {
                fn default() -> Self {
                    Self {
                        $(
                            $in: Default::default(),
                        )*
                        $($param_name: Default::default()),*
                    }
                }
            }
            impl Circuit for $name {
                fn exec<B: Backend>(&self, frontend: &Frontend<B>) {
                    let $executor = frontend;
                    $(
                        let $in = self.$in;
                    )*
                    $(
                        let $param_name = self.$param_name;
                    )*
                    $(
                        let $out;
                    )*
                    $body
                    $(
                        $executor.output($out);
                    )*
                }
            }
        }
    };
}

#[doc(inline)]
pub use _define_test_circuit as define_test_circuit;

#[macro_export]
#[doc(hidden)]
macro_rules! _define_test_func {
    (
        $func: ident,
        $num_samples: expr,
        {$($in:ident : $in_t:ty),* $(,)?},
        [$($out:ident),* $(,)?],
        |$executor: ident| $body: block,
        $ref_body: block
        $(, $param_name: ident : $param_type: ty, $param_values: expr)?
    ) => {
        ::paste::paste! {
            fn $func(){
                use zkboo::{
                    circuit::Circuit,
                    backend::{Backend, Frontend},
                    executor::{exec, OwnedFlexibleWordPool},
                    word::{Words, WordLike},
                };
                use $crate::common::{rand_words::test_vec, test_func::define_test_circuit, proofs::test_proof};
                type WP = OwnedFlexibleWordPool<usize>;
                define_test_circuit!(TestCircuit, {$($in: $in_t,)*}, [$($out,)*], |$executor| $body $(, $param_name: $param_type)*);
                let _param_values = [1usize; 1];
                $(
                    let _param_values = $param_values;
                )*
                for _p in _param_values {
                    $(
                        let $param_name = _p;
                    )*
                    let mut seed = 0u64;
                    $(
                        seed += 1;
                        let mut [<iter_ $in>] = test_vec::<_, _, $in_t>($num_samples, seed).into_iter();
                    )*
                    for _  in 0..$num_samples {
                        $(
                            let [<_ $in>] = [<iter_ $in>].next().unwrap();
                        )*
                        // Test execution:
                        let circuit = TestCircuit {$($in: [<_ $in>],)* $($param_name,)*};
                        let outputs = exec::<_, WP>(&circuit);
                        // Reference execution:
                        $(
                            let $in = [<_ $in>];
                        )*
                        $(
                            let $out;
                        )*
                        $ref_body
                        let mut expected_outputs = Words::new();
                        $(
                            expected_outputs.as_vec_mut().extend(WordLike::to_word($out).to_le_words());
                        )*
                        // Ensure matching output words:
                        assert_eq!(outputs, expected_outputs);
                        // Test proof generation:
                        test_proof(&circuit);
                    }
                }
            }
        }
    };
}

#[doc(inline)]
pub use _define_test_func as define_test_func;

#[macro_export]
#[doc(hidden)]
macro_rules! _run_test {
    (
        $func: ident,
        $num_samples: expr,
        {$($in:ident : $in_t:ty),* $(,)?},
        [$($out:ident),* $(,)?],
        |$executor: ident| $body: block,
        $ref_body: block
        $(, $param_name: ident : $param_type: ty, $param_values: expr)?
    ) => {
        use $crate::common::test_func::define_test_func;
        define_test_func!(
            $func,
            $num_samples,
            {$($in : $in_t),*},
            [$($out),*],
            |$executor| $body,
            $ref_body
            $(, $param_name : $param_type, $param_values)*
        );
        $func();
    };
}

#[doc(inline)]
pub use _run_test as run_test;
