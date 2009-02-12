#include <utility>
template<typename TYPE_FLOAT>
boost::numeric::ublas::matrix<TYPE_FLOAT>    quaternion_to_R3_rotation(::boost::math::quaternion<TYPE_FLOAT> const & q)
{
    using    ::std::numeric_limits;
    
    TYPE_FLOAT    a = q.R_component_1();
    TYPE_FLOAT    b = q.R_component_2();
    TYPE_FLOAT    c = q.R_component_3();
    TYPE_FLOAT    d = q.R_component_4();
    
    TYPE_FLOAT    aa = a*a;
    TYPE_FLOAT    ab = a*b;
    TYPE_FLOAT    ac = a*c;
    TYPE_FLOAT    ad = a*d;
    TYPE_FLOAT    bb = b*b;
    TYPE_FLOAT    bc = b*c;
    TYPE_FLOAT    bd = b*d;
    TYPE_FLOAT    cc = c*c;
    TYPE_FLOAT    cd = c*d;
    TYPE_FLOAT    dd = d*d;
    
    TYPE_FLOAT    norme_carre = aa+bb+cc+dd;
    
    if    (norme_carre <= numeric_limits<TYPE_FLOAT>::epsilon())
    {
        ::std::string            error_reporting("Argument to quaternion_to_R3_rotation is too small!");
        ::std::underflow_error   bad_argument(error_reporting);
        
        throw(bad_argument);
    }
    
    boost::numeric::ublas::matrix<TYPE_FLOAT>    out_matrix(3,3);
    
    out_matrix(0,0) = (aa+bb-cc-dd)/norme_carre;
    out_matrix(0,1) = 2*(-ad+bc)/norme_carre;
    out_matrix(0,2) = 2*(ac+bd)/norme_carre;
    out_matrix(1,0) = 2*(ad+bc)/norme_carre;
    out_matrix(1,1) = (aa-bb+cc-dd)/norme_carre;
    out_matrix(1,2) = 2*(-ab+cd)/norme_carre;
    out_matrix(2,0) = 2*(-ac+bd)/norme_carre;
    out_matrix(2,1) = 2*(ab+cd)/norme_carre;
    out_matrix(2,2) = (aa-bb-cc+dd)/norme_carre;
    
    return(out_matrix);
}

