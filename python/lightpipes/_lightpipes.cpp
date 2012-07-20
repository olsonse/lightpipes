
#include <lightpipes/Field.h>

#include <boost/python.hpp>
#include <numpy/noprefix.h>

#include <ostream>
#include <sstream>
#include <iostream>

using namespace boost::python;

namespace lightpipes {
  std::ostream & operator<< ( std::ostream & out, const Field::Info & i ) {
    out << "Info(number=" << i.number << ", "
             "side_length= " << i.side_length << ", "
             "wavelength= " << i.lambda << ", "
             "fft_level= " << i.fft_level << ", "
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
    npy_intp dims[2] = {f.info.number, f.info.number};
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
                      const std::string & X_filename = "",
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

    f.steps( step_size, number_steps, n_data, X_filename, dump_period );
  }

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    Steps, steps_wrapper, 2, 6)

  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    Interpolate, Field::interpolate, 0, 6)


  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    CircularAperture, Field::circular_aperture, 1, 3)

  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    CircularScreen, Field::circular_screen, 1, 3)


  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    RectangularAperture, Field::rectangular_aperture, 1, 5)

  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    RectangularScreen, Field::rectangular_screen, 1, 5)


  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    GaussianAperture, Field::gaussian_aperture, 1, 4)

  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    GaussianScreen, Field::gaussian_screen, 1, 4)


  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    SuperGaussianAperture, Field::supergaussian_aperture, 2, 5)

  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    SuperGaussianScreen, Field::supergaussian_screen, 2, 5)

}// namespace (anonymous)

BOOST_PYTHON_MODULE(_lightpipes) {
  import_array(); // for importing numpy array stuff
  numeric::array::set_module_and_type("numpy", "ndarray");

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
    init<unsigned int, double, double,
         optional< std::complex<double>, int, double > >(
      args("number","side_length","lambda",
           "init","fft_level","sph_coords_factor"),
      " number      : Number of side-elements stored in field.      \n"
      "               Note:  sizeof(Field) = number^2.              \n"
      " side_length : Physical size of the side of the field.       \n"
      " lambda      : Wavelength of the field.                      \n"
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
                                                            return_self<>())
    .def("lens",                        &Field::lens,       return_self<>() )
    .def("t_lens",                      &Field::t_lens,     return_self<>() )
    .def("axicon",                      &Field::axicon,     return_self<>() )
    .def("fresnel",                     &Field::fresnel,    return_self<>() )
    .def("forward",                     &Field::forward,    return_self<>() )
    .def("lens_forvard",                &Field::lens_forvard, return_self<>() )
    .def("lens_fresnel",                &Field::lens_fresnel, return_self<>() )
    .def("forvard",                     &Field::forvard,    return_self<>() )
    .def("spherical_to_normal_coords",  &Field::spherical_to_normal_coords,
                                                            return_self<>() )
    .def("circular_aperture",           &Field::circular_aperture,
      CircularAperture( args("r","x0","y0"),
                        "Defaults:  x0=0.0, y0=0.0" )     [ return_self<>() ] )
    .def("circular_screen",             &Field::circular_screen,
      CircularScreen( args("r","x0","y0"),
                      "Defaults:  x0=0.0, y0=0.0" )       [ return_self<>() ] )
    .def("rectangular_aperture",        &Field::rectangular_aperture,
      RectangularAperture( args("Lx","Ly","x0","y0","angle"),
        "Defaults: Ly=-0.1, x0=0.0, y0=0.0, angle=0.0" )  [ return_self<>() ] )
    .def("rectangular_screen",          &Field::rectangular_screen,
      RectangularScreen( args("Lx","Ly","x0","y0","angle"),
        "Defaults: Ly=-0.1, x0=0.0, y0=0.0, angle=0.0" )  [ return_self<>() ] )
    .def("supergaussian_aperture",      &Field::supergaussian_aperture,
      SuperGaussianAperture( args("w","n","x0","y0","A"),
        "Defaults: x0=0.0, y0=0.0, A=1.0" )               [ return_self<>() ] )
    .def("supergaussian_screen",        &Field::supergaussian_screen,
      SuperGaussianScreen( args("w","n","x0","y0","A"),
        "Defaults: x0=0.0, y0=0.0, A=1.0" )               [ return_self<>() ] )
    .def("gaussian_aperture",           &Field::gaussian_aperture,
      GaussianAperture( args("w","x0","y0","A"),
        "Defaults: x0=0.0, y0=0.0, A=1.0" )               [ return_self<>() ] )
    .def("gaussian_screen",             &Field::gaussian_screen,
      GaussianScreen( args("w","x0","y0","A"),
        "Defaults: x0=0.0, y0=0.0, A=1.0" )               [ return_self<>() ] )
    .def("fft3",                        &Field::fft3,       return_self<>() )
    .def("tilt",                        &Field::tilt,       return_self<>() )
    .def("zernike",                     &Field::zernike,    return_self<>() )
    .def("get_strehl",                  &Field::get_strehl                  )
    .def("pip_fft",                     &Field::pip_fft,    return_self<>() )
    .def("compatible",                  &Field::compatible                  )
    .def("normalize",                   &Field::normalize,  return_self<>() )
    .def("l_amplify",                   &Field::l_amplify,  return_self<>() )
    .def("interpolate",                 &Field::interpolate,
      Interpolate( args("new_side_length", "new_number",
                        "x_shift", "y_shift", "angle", "magnif"),
        "Defaults: new_side_length=0.0\n"
        "          new_number     =0  \n"
        "          x_shift        =0.0\n"
        "          y_shift        =0.0\n"
        "          angle          =0.0\n"
        "          magnif         =1.0"
      )[ return_self<>() ] )
    .def("steps",                       steps_wrapper,
      Steps(args("step_size","N","n","X_filename","dump_period"),
        " step_size       \n"
        "   Propagation distance to make for a single step\n"
        " N          [=1] \n"
        "   Number of steps of step_size to make\n"
        " n [=ones(info.number,info.number) * (1+0j)]\n"
        "   Complex index of refraction as a function of position\n"
        " X_filename [='']\n"
        " dump_period[=1]"
      )[ return_self<>() ]
    )
    .def("copy",                        &Field::copy)
    ;
}

