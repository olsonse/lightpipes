
#include <lightpipes/Field.h>

#include <boost/python.hpp>
#include <numpy/noprefix.h>

#include <ostream>
#include <sstream>
#include <iostream>

#include <cassert>

using namespace boost::python;

namespace lightpipes {
  std::ostream & operator<< ( std::ostream & out, const Field::Info & i ) {
    out << "Info(number=" << i.number << ", "
             "side_length=" << i.side_length << ", "
             "wavelength=" << i.lambda << ", "
             "fft_level=" << i.fft_level << ", "
             "sph_coords_factor= " << i.sph_coords_factor << ")";
    return out;
  }
}

namespace { // (anonymous) namespace

  template <typename T>
  std::string to_string( const T & t ) {
    std::ostringstream ostr;
    ostr << t;
    return ostr.str();
  }

  using namespace lightpipes;


  PyObject * make_ndarray_wrapper( Field & f ) {
    assert( f.info.number.first <= (2uL<<(8*sizeof(f.info.number.first)-1))-1 );
    assert( f.info.number.second <= (2uL<<(8*sizeof(f.info.number.second)-1))-1 );
    npy_long dims[2] = {(long)f.info.number.first, (long)f.info.number.second};
    return PyArray_SimpleNewFromData(2, dims, PyArray_CDOUBLE, f.val);
  }

  numeric::array empty_array() {
    int zero = 0;
    return
      extract< numeric::array >(
        object( handle<>( PyArray_FromDims(0, &zero, PyArray_CDOUBLE) ) )
      );
  }

  void steps_wrapper( Field & f,
                      const double & step_size,
                      const int & number_steps = 1,
                      const numeric::array & n = empty_array(),
                      const std::string & dump_filename = "",
                      const int & dump_period = 1 ) {
    typedef numeric::array A;
    typedef std::complex<double> Complex;
    Complex * n_data = NULL;
    A tmp( empty_array() ); //keep in scope until end of function

    if ( ! PyArray_IsZeroDim( n.ptr() ) ) {
      if ( PyArray_TYPE( n.ptr() ) != PyArray_CDOUBLE ) {
        // PyArray_CDOUBLE == 'D'
        tmp = static_cast<A>(static_cast<A>(n.copy()).astype('D'));
        n_data = static_cast<Complex*>(PyArray_DATA( tmp.ptr() ));
      } else {
        n_data = static_cast<Complex*>(PyArray_DATA(n.ptr()));
      }
    }

    f.steps( step_size, number_steps, n_data, dump_filename, dump_period );
  }

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    Steps, steps_wrapper, 2, 6)


  template<class T1, class T2=T1>
  struct PairToTupleConverter {
    static PyObject* convert(const Pair<T1, T2>& pair) {
      return incref(make_tuple(pair.first, pair.second).ptr());
    }
  };

  template<class T1, class T2=T1>
  struct StdPairToTupleConverter {
    static PyObject* convert(const std::pair<T1, T2>& pair) {
      return incref(make_tuple(pair.first, pair.second).ptr());
    }
  };

}// namespace (anonymous)

#define class_Pair(type,name) \
  class_< Pair<type> >(name) \
    .def(init< type, optional< type > >( \
           args("first","second"), \
           " first : first component of pair \n" \
           " second: second component of pair \n" \
           "         [Default : <first>]\n" )) \
    .def_readonly("first",  &Pair<type>::first) \
    .def_readonly("second", &Pair<type>::second) \
    .def( "__str__", to_string< Pair<type> > ) \
    .def( "__repr__", to_string< Pair<type> > )

BOOST_PYTHON_MODULE(_lightpipes) {
  import_array(); // for importing numpy array stuff
  numeric::array::set_module_and_type("numpy", "ndarray");
  scope().attr("version") = LIGHTPIPES_VERSION;

  //to_python_converter<Pair<size_t>, PairToTupleConverter<size_t> >();
  //to_python_converter<Pair<double>, PairToTupleConverter<double> >();
  //to_python_converter<std::pair<size_t,size_t>,
  //  StdPairToTupleConverter<size_t,size_t> >();
  //to_python_converter<std::pair<double,double>,
  //  StdPairToTupleConverter<double,double> >();

  class_Pair(int,"IntPair");
  class_Pair(unsigned int,"UIntPair");
  class_Pair(long,"LongPair");
  class_Pair(size_t,"SizePair");
  class_Pair(double,"DoublePair");


  class_<Field::Info>("Info")
    .def_readonly("number", &Field::Info::number)
    .def_readwrite("side_length", &Field::Info::side_length)
    .def_readwrite("wavelength", &Field::Info::lambda)
    .def_readonly("fft_level", &Field::Info::fft_level)
    .def_readonly("sph_coords_factor", &Field::Info::sph_coords_factor)
    .def("compatible", &Field::Info::compatible)
    .def( "__str__", to_string<Field::Info> )
    .def( "__repr__", to_string<Field::Info> )
    ;

  class_<Field>(
    "Field",
    init<Pair<size_t>, Pair<double>, double,
         optional< std::complex<double>, int, double > >(
      args("number","side_length","wavelength",
           "init","fft_level","sph_coords_factor"),
      " number      : Number of side-elements stored in field.      \n"
      "               Note:  sizeof(Field) = number^2.              \n"
      " side_length : Physical size of the side of the field.       \n"
      " wavelength  : Wavelength of the field.                      \n"
      " init        : Initial (constant) value of field             \n"
      "               [Default : 1.0 + 0j]                          \n"
      " fft_level   : FFT Status [Default : 0].                     \n"
      " sph_coords_factor : Spherical coordinates conversion factor.\n"
      "               [Default : 0.0]                               \n"
    )
  )
    .def(init<Field::Info, optional< std::complex<double> > >(
           args("info","init"),
           " info        : Initialize based on Field.info.               \n"
           " init        : Initial (constant) value of field             \n"
           "               [Default : 1.0 + 0j]                          \n"
         ))
    //.def_readonly("number", &Field::Info::number),
    //.def_readwrite("side_length", &Field::Info::side_length),
    .add_property("info", 
      make_function( &Field::getinfo, return_internal_reference<>()) )
    .add_property("value", make_ndarray_wrapper)
    .def(self + self)
    .def(self += self)
    .def(self *= float())
    .def(self *= int())
    .def(self *= double())
    .def(self *= std::complex<double>())
    .def("fill",                        &Field::fill<double>, return_self<>())
    .def("fill",                        &Field::fill< std::complex<double> >,
                                                            return_self<>() )
    .def("lens",                        &Field::lens,
         ( arg("f"), arg("x0")=0.0, arg("y0")=0.0 ),        return_self<>() )
    .def("t_lens",                      &Field::t_lens,
         (arg("R"),arg("f"),arg("x0")=0.0,arg("y0")=0.0),   return_self<>() )
    .def("axicon",                      &Field::axicon,
         (arg("r"),arg("n1")=std::complex<double>(1.5,0.),
          arg("x0")=0.0,arg("y0")=0.0 ),                    return_self<>() )
    .def("fresnel",                     &Field::fresnel,
         (arg("z")),                                        return_self<>(),
          "Propagates Field using convolution.\n"
          " z\n"
          "   Distance to propagate.\n")
    .def("forward",                     &Field::forward,
         (arg("z"),arg("new_side_length"),arg("new_number")),return_self<>(),
         "Propagates Field using direct integration.\n")
    .def("lens_forvard",                &Field::lens_forvard,
         (arg("f"),arg("z")),                               return_self<>(),
         "Propagates Field in spherical coordinates using FFT.            \n"
         " f                        \n"
         "   Focal distance of lens that determines the curvature of the  \n"
         "   coordinate system.     \n"
         " z                        \n"
         "   Propagation distance.  \n")
    .def("lens_fresnel",                &Field::lens_fresnel,
         (arg("f"),arg("z")),                               return_self<>(),
         "Propagates Field in spherical coordinates.                      \n"
         " f                        \n"
         "   Focal distance of lens that determines the curvature of the  \n"
         "   coordinate system.     \n"
         " z                        \n"
         "   Propagation distance.  \n")
    .def("forvard",                     &Field::forvard,
         (arg("z")),                                        return_self<>(),
          "Propagates Field using FFT. \n"
          " z\n"
          "   Distance to propagate.\n")
    .def("spherical_to_normal_coords",  &Field::spherical_to_normal_coords,
                                                            return_self<>() )
    .def("circular_aperture",           &Field::circular_aperture,
         (arg("r"),arg("x0")=0.0,arg("y0")=0.0),            return_self<>() )
    .def("circular_screen",             &Field::circular_screen,
         (arg("r"),arg("x0")=0.0,arg("y0")=0.0),            return_self<>() )
    .def("rectangular_aperture",        &Field::rectangular_aperture,
         (arg("Lx"),arg("Ly")=-0.1,
          arg("x0")=0.0,arg("y0")=0.0,arg("angle")=0.0),    return_self<>() )
    .def("rectangular_screen",          &Field::rectangular_screen,
         (arg("Lx"),arg("Ly")=-0.1,
          arg("x0")=0.0,arg("y0")=0.0,arg("angle")=0.0),    return_self<>() )
    .def("supergaussian_aperture",      &Field::supergaussian_aperture,
         (arg("wx"),arg("wy"),arg("n"),
          arg("x0")=0.0,arg("y0")=0.0,arg("A")=1.0),        return_self<>() )
    .def("supergaussian_screen",        &Field::supergaussian_screen,
         (arg("wx"),arg("wy"),arg("n"),
          arg("x0")=0.0, arg("y0")=0.0, arg("A")=1.0),      return_self<>() )
    .def("gaussian_aperture",           &Field::gaussian_aperture,
         (arg("wx"),arg("wy")=0.0,arg("x0")=0.0,arg("y0")=0.0,arg("A")=1.0),
          return_self<>() )
    .def("gaussian_screen",             &Field::gaussian_screen,
         (arg("wx"),arg("wy")=0.0,arg("x0")=0.0,arg("y0")=0.0,arg("A")=1.0),
          return_self<>() )
    .def("fft3",                        &Field::fft3,       return_self<>() )
    .def("tilt",                        &Field::tilt,       return_self<>() )
    .def("zernike",                     &Field::zernike,
         (arg("n"),arg("m"),arg("R"),arg("A")),             return_self<>(),
          "Zernike polynomials.\n"
          " n\n"
          "   Radial order\n"
          " m\n"
          "   Azimuthal order\n"
          " R\n"
          "   Radius of amplitute spec\n"
          " A\n"
          "   Amplitude of abberation at r=R, in radians\n")
    .def("get_strehl",                  &Field::get_strehl                  )
    .def("pip_fft",                     &Field::pip_fft,    return_self<>() )
    .def("compatible",                  &Field::compatible                  )
    .def("normalize",                   &Field::normalize,  return_self<>() )
    .def("l_amplify",                   &Field::l_amplify,  return_self<>() )
    .def("interpolate",                 &Field::interpolate,
         (arg("new_side_length")=0.0, arg("new_number")=0,
          arg("x_shift")=0.0, arg("y_shift")=0.0,
          arg("angle")=0.0, arg("magnif")=1.0),             return_self<>() )
    .def("steps",                       steps_wrapper,
      Steps(args("step_size","N","n","dump_filename","dump_period"),
        " Propagate the field using finite-difference routine.\n"
        "\n"
        " step_size       \n"
        "   Propagation distance to make for a single step\n"
        " N          [=1] \n"
        "   Number of steps of step_size to make\n"
        " n [=ones(info.number,info.number) * (1+0j)]\n"
        "   Complex index of refraction as a function of position\n"
        " dump_filename [='']\n"
        " dump_period[=1]"
      )[ return_self<>() ]
    )
    .def("copy",                        &Field::copy)
    ;
}

