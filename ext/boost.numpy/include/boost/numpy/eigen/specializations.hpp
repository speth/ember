#ifndef BOOST_PYTHON_EIGEN_SPECIALIZATIONS_HPP_INCLUDED
#define BOOST_PYTHON_EIGEN_SPECIALIZATIONS_HPP_INCLUDED

#include <boost/numpy/eigen/NumpyView.hpp>
#include <boost/numpy/eigen/eigen_to_python.hpp>
#include <boost/numpy/eigen/return_internal_matrix.hpp>

#define BOOST_PYTHON_EIGEN_TO_PYTHON_VALUE(ARGS, TYPE, CV)              \
    namespace boost { namespace python {                                \
    template < ARGS >                                                   \
    struct to_python_value< TYPE CV > : public detail::builtin_to_python { \
        inline PyObject * operator()(TYPE CV x) const {                 \
            try {                                                       \
                object r = eigen_to_python< TYPE >::to_python_value(x); \
                Py_INCREF(r.ptr());                                     \
                return r.ptr();                                         \
            } catch (error_already_set & exc) {                         \
                handle_exception();                                     \
                return NULL;                                            \
            }                                                           \
        }                                                               \
        inline PyTypeObject const * get_pytype() const {                \
            return converter::object_manager_traits<numpy::matrix>::get_pytype(); \
        }                                                               \
    };                                                                  \
    }}

#define BOOST_PYTHON_EIGEN_ARG_TO_PYTHON(ARGS, TYPE)                    \
    namespace boost { namespace python { namespace converter {          \
        template < ARGS >                                               \
        struct arg_to_python< TYPE > : public handle<> {                \
            inline arg_to_python(TYPE const & v) :                      \
                handle<>(python::to_python_value<TYPE const &>()(v)) {} \
        };                                                              \
    }}}

#define BOOST_PYTHON_EIGEN_FROM_PYTHON_AUTO(ARGS, TYPE)        \
    namespace boost { namespace python { namespace converter {          \
    template < ARGS >                                                   \
    struct arg_rvalue_from_python< TYPE const & > {                     \
        typedef TYPE result_type;                                       \
        typedef typename TYPE::Scalar _Scalar;                          \
        typedef Eigen::NumpyView<_Scalar, TYPE::RowsAtCompileTime,TYPE::ColsAtCompileTime, \
                                 typename Eigen::internal::traits< TYPE >::XprKind> NumpyViewType; \
        arg_rvalue_from_python(PyObject * p) : _p(borrowed(p)) {}       \
        bool convertible() const {                                      \
            try {                                                       \
                object arg(_p);                                         \
                numpy::matrix m(arg, numpy::dtype::get_builtin<_Scalar>()); \
                NumpyViewType pm(m);                                          \
                _p = handle<>(borrowed(m.ptr()));                       \
            } catch (error_already_set & exc) {                         \
                PyErr_Clear();                                          \
                return false;                                           \
            }                                                           \
            return true;                                                \
        }                                                               \
        result_type operator()() const {                                \
            boost::numpy::matrix m(python::detail::borrowed_reference(_p.get())); \
            NumpyViewType pm(m);                                        \
            return result_type(pm);                                     \
        }                                                               \
    private:                                                            \
        mutable handle<> _p;                                            \
    };                                                                  \
    }}}

#define EIGEN_MATRIX_ARGS typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols
#define EIGEN_MATRIX_TYPE Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>
BOOST_PYTHON_EIGEN_TO_PYTHON_VALUE(EIGEN_MATRIX_ARGS, EIGEN_MATRIX_TYPE, &)
BOOST_PYTHON_EIGEN_TO_PYTHON_VALUE(EIGEN_MATRIX_ARGS, EIGEN_MATRIX_TYPE, const &)
BOOST_PYTHON_EIGEN_ARG_TO_PYTHON(EIGEN_MATRIX_ARGS, EIGEN_MATRIX_TYPE)
BOOST_PYTHON_EIGEN_FROM_PYTHON_AUTO(EIGEN_MATRIX_ARGS, EIGEN_MATRIX_TYPE)
#undef EIGEN_MATRIX_ARGS
#undef EIGEN_MATRIX_TYPE

#define EIGEN_ARRAY_ARGS typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols
#define EIGEN_ARRAY_TYPE Eigen::Array<Scalar,Rows,Cols,Options,MaxRows,MaxCols>
BOOST_PYTHON_EIGEN_TO_PYTHON_VALUE(EIGEN_ARRAY_ARGS, EIGEN_ARRAY_TYPE, &)
BOOST_PYTHON_EIGEN_TO_PYTHON_VALUE(EIGEN_ARRAY_ARGS, EIGEN_ARRAY_TYPE, const &)
BOOST_PYTHON_EIGEN_ARG_TO_PYTHON(EIGEN_ARRAY_ARGS, EIGEN_ARRAY_TYPE)
BOOST_PYTHON_EIGEN_FROM_PYTHON_AUTO(EIGEN_ARRAY_ARGS, EIGEN_ARRAY_TYPE)
#undef EIGEN_ARRAY_ARGS
#undef EIGEN_ARRAY_TYPE

#endif // !BOOST_PYTHON_EIGEN_SPECIALIZATIONS_HPP_INCLUDED
