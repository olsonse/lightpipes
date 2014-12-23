
#ifndef OCT_FIELD_H
#define OCT_FIELD_H

#include <octave/oct.h>
#include <octave/ov-struct.h>
#include <lightpipes/Field.h>

using lightpipes::Field;

Field::Info & mapToInfo(const Octave_map & ov_info, Field::Info & info) {
    info.number.first           = ov_info.contents("number_x")(0).int_value();
    info.number.second          = ov_info.contents("number_y")(0).int_value();
    info.side_length.first      = ov_info.contents("side_length_x")(0).double_value();;
    info.side_length.second     = ov_info.contents("side_length_y")(0).double_value();;
    info.lambda                 = ov_info.contents("lambda")(0).double_value();;
    info.fft_level              = ov_info.contents("fft_level")(0).int_value();;
    info.sph_coords_factor      = ov_info.contents("sph_coords_factor")(0).double_value();;
    return info;
}

Octave_map infoToMap(const Field::Info & info) {
    Octave_map ov_info;
    ov_info.assign("number_x", info.number.first);
    ov_info.assign("number_y", info.number.second);
    ov_info.assign("side_length_x", info.side_length.first);
    ov_info.assign("side_length_y", info.side_length.second);
    ov_info.assign("lambda", info.lambda);
    ov_info.assign("fft_level", info.fft_level);
    ov_info.assign("sph_coords_factor", info.sph_coords_factor);
    return ov_info;
}

bool mapIsValidInfo( std::ostream & errout, const Octave_map & ov_info ) {
    if ( !ov_info.contains("number_x") ||
         !ov_info.contains("number_y") ||
         !ov_info.contains("side_length_x") ||
         !ov_info.contains("side_length_y") ||
         !ov_info.contains("lambda") ||
         !ov_info.contains("fft_level") ||
         !ov_info.contains("sph_coords_factor") ) {
        errout << "Field info invalid" << std::endl;
        return false;
    } else return true;
}

bool mapIsValidField (std::ostream & errout, const Octave_map & Fmap) {
    if ( !Fmap.contains("info") ||
         !mapIsValidInfo(errout, Fmap.contents("info")(0).map_value()) ||
         !Fmap.contains("F") ||
         !Fmap.contents("F")(0).is_matrix_type() ) {
        errout << "Field invalid" << std::endl;
        return false;
    } else return true;
}

Field * mapToField( const Octave_map  & Fmap ) {
    Field::Info info;
    Octave_map ov_info = Fmap.contents("info")(0).map_value();
    ComplexMatrix F = Fmap.contents("F")(0).complex_matrix_value();
    Field * field = new Field(mapToInfo(ov_info, info));

    int ij = 0;
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < F.columns(); j++) {
            field->val[ij++] = F(i,j);
        }/* for */
    }/* for */

    return field;
}

ComplexMatrix fieldToCMatrix( const Field & field ) {
    ComplexMatrix F(field.info.number.first, field.info.number.second);

    int ij = 0;
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < F.columns(); j++) {
            F(i,j) = field[ij++];
        }/* for */
    }/* for */

    return F;
}

Octave_map fieldToMap ( const Field & field ) {
    Octave_map retval;
    retval.assign("F", fieldToCMatrix(field));
    retval.assign("info", infoToMap(field.info));

    return retval;
}

#endif // OCT_FIELD_H
