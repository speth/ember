#ifndef BOOST_PYTHON_EIGEN_NUMPYVIEW_HPP_INCLUDED
#define BOOST_PYTHON_EIGEN_NUMPYVIEW_HPP_INCLUDED

#include <boost/numpy.hpp>
#include <Eigen/Core>

namespace Eigen {

template <typename _Scalar, int _Rows, int _Cols, typename _XprKind> class NumpyView;

namespace internal {

template <typename _Scalar, int _Rows, int _Cols, typename _XprKind>
struct traits< NumpyView<_Scalar,_Rows,_Cols,_XprKind> > {
    typedef _Scalar Scalar;
    typedef int Index;
    typedef Dense StorageKind;
    typedef _XprKind XprKind;
    enum {
        RowsAtCompileTime = _Rows,
        ColsAtCompileTime = _Cols,
        MaxRowsAtCompileTime = _Rows,
        MaxColsAtCompileTime = _Cols,
        InnerStrideAtCompileTime = Eigen::Dynamic,
        OuterStrideAtCompileTime = Eigen::Dynamic,
        IsVectorAtCompileTime = 1,
        Flags = Eigen::DirectAccessBit | Eigen::NestByRefBit | Eigen::RowMajorBit | Eigen::LvalueBit,
        CoeffReadCost = NumTraits<Scalar>::ReadCost
    };
};

} // namespace internal

template <typename _Scalar, int _Rows, int _Cols, typename _XprKind>
class NumpyView : public internal::dense_xpr_base< NumpyView<_Scalar,_Rows,_Cols,_XprKind> >::type {
public:

    typedef typename internal::dense_xpr_base< NumpyView<_Scalar,_Rows,_Cols,_XprKind> >::type Base;

    EIGEN_GENERIC_PUBLIC_INTERFACE(NumpyView);

    EIGEN_INHERIT_ASSIGNMENT_OPERATORS(NumpyView);

    explicit NumpyView(boost::numpy::matrix const & py) : _py(py) {
        if (
            (_Rows != Eigen::Dynamic && _Rows != _py.shape(0))
            || (_Cols != Eigen::Dynamic && _Cols != _py.shape(1))
        ) {
            if (
                (_Rows == 1 && _py.shape(1) == 1 && (_Cols == Eigen::Dynamic || _Cols == _py.shape(0)))
                || (_Cols == 1 && _py.shape(0) == 1 && (_Rows == Eigen::Dynamic || _Rows == _py.shape(1)))
            ) {
                _py = py.transpose();
            } else {
                PyErr_SetString(PyExc_ValueError, "Incorrect shape for fixed-size matrix.");
                boost::python::throw_error_already_set();
            }
        }
        if (_py.get_dtype() != boost::numpy::dtype::get_builtin<_Scalar>()) {
            PyErr_SetString(PyExc_TypeError, "Incorrect data type for matrix.");
            boost::python::throw_error_already_set();
        }
    }

    inline Index innerStride() const { return _py.strides(1); }
    inline Index outerStride() const { return _py.strides(0); }
    inline Index rowStride() const { return _py.strides(0); }
    inline Index colStride() const { return _py.strides(1); }

    inline int rows() const { return _py.shape(0); }
    inline int cols() const { return _py.shape(1); }

    Scalar * data() { return _py.get_data(); }
    Scalar const * data() const { return _py.get_data(); }

    inline Scalar & coeffRef(int row, int col) {
        return *reinterpret_cast<Scalar*>(_py.get_data() + row * _py.strides(0) + col * _py.strides(1));
    }

    inline Scalar & coeffRef(int index) {
        if (_Rows == 1)
            return *reinterpret_cast<Scalar*>(_py.get_data() + index * _py.strides(1));
        else if (_Cols == 1)
            return *reinterpret_cast<Scalar*>(_py.get_data() + index * _py.strides(0));
        PyErr_SetString(PyExc_TypeError, "Cannot index NumpyView with a single integer.");
        boost::python::throw_error_already_set();
    }

    inline Scalar const & coeffRef(int row, int col) const {
        return *reinterpret_cast<Scalar*>(_py.get_data() + row * _py.strides(0) + col * _py.strides(1));
    }

    inline Scalar const & coeffRef(int index) const {
        if (_Rows == 1)
            return *reinterpret_cast<Scalar*>(_py.get_data() + index * _py.strides(1));
        else if (_Cols == 1)
            return *reinterpret_cast<Scalar*>(_py.get_data() + index * _py.strides(0));
        PyErr_SetString(PyExc_TypeError, "Cannot index NumpyView with a single integer.");
        boost::python::throw_error_already_set();
    }

    inline Scalar const & coeff(int row, int col) const {
        return *reinterpret_cast<Scalar*>(_py.get_data() + row * _py.strides(0) + col * _py.strides(1));
    }

    inline Scalar const & coeff(int index) const {
        if (_Rows == 1)
            return *reinterpret_cast<Scalar*>(_py.get_data() + index * _py.strides(1));
        else if (_Cols == 1)
            return *reinterpret_cast<Scalar*>(_py.get_data() + index * _py.strides(0));
        PyErr_SetString(PyExc_TypeError, "Cannot index NumpyView with a single integer.");
        boost::python::throw_error_already_set();
    }

    boost::numpy::matrix getPyObject() const { return _py; }

private:
    boost::numpy::matrix _py;
};

} // namespace Eigen

namespace boost { namespace python {

template <typename _Scalar, int _Rows, int _Cols, typename _XprKind>
struct to_python_value< Eigen::NumpyView<_Scalar,_Rows,_Cols,_XprKind> const & >
    : public detail::builtin_to_python
{
    inline PyObject * operator()(Eigen::NumpyView<_Scalar,_Rows,_Cols,_XprKind> const & x) const {
        numpy::matrix obj(x.getPyObject());
        Py_INCREF(obj.ptr());
        return obj.ptr();
    }
    inline PyTypeObject const * get_pytype() const {
        return converter::object_manager_traits<numpy::matrix>::get_pytype();
    }
};

template <typename _Scalar, int _Rows, int _Cols, typename _XprKind>
struct to_python_value< Eigen::NumpyView<_Scalar,_Rows,_Cols,_XprKind> & >
    : public detail::builtin_to_python
{
    inline PyObject * operator()(Eigen::NumpyView<_Scalar,_Rows,_Cols,_XprKind> & x) const {
        numpy::matrix obj(x.getPyObject());
        Py_INCREF(obj.ptr());
        return obj.ptr();
    }
    inline PyTypeObject const * get_pytype() const {
        return converter::object_manager_traits<numpy::matrix>::get_pytype();
    }
};

namespace converter {

template <typename _Scalar, int _Rows, int _Cols, typename _XprKind>
struct arg_to_python< Eigen::NumpyView<_Scalar,_Rows,_Cols,_XprKind> > : public handle<> {
    inline arg_to_python(Eigen::NumpyView<_Scalar,_Rows,_Cols,_XprKind> const & v) :
        handle<>(borrowed(v.getPyObject().ptr())) {}
};

template <typename _Scalar, int _Rows, int _Cols, typename _XprKind>
struct arg_rvalue_from_python< Eigen::NumpyView<_Scalar,_Rows,_Cols,_XprKind> const & > {
    typedef Eigen::NumpyView<_Scalar,_Rows,_Cols,_XprKind> result_type;

    arg_rvalue_from_python(PyObject * p) : _p(borrowed(p)) {}

    bool convertible() const { return true; }

    result_type operator()() const {
        boost::numpy::matrix m(python::detail::borrowed_reference(_p.get()));
        return Eigen::NumpyView<_Scalar,_Rows,_Cols,_XprKind>(m);
    }

private:
    mutable handle<> _p;
};

} // namespace boost::python::converter

}} // namespace boost::python

#endif // !BOOST_PYTHON_EIGEN_NUMPYVIEW_HPP_INCLUDED
