%DIMENSIONS.M - Script used by STRUCT2D.M
  re_name{1} = 'S1';
    re_xc{1} = '-sep';
    re_yc{1} = 'd+th1/2';
 re_shape{1} = 'rectangle';
 re_width{1} = 'w';
re_height{1} = 'th1';
re_diel(1) = false;
re_diel_r{1} = '';

  re_name{2} = 'S2';
    re_xc{2} = 'sep';
    re_yc{2} = 'd+th1/2';
 re_shape{2} = 'rectangle';
 re_width{2} = 'w';
re_height{2} = 'th1';
re_diel(2) = false;
re_diel_r{2} = '';

  re_name{3} = 'S3';
    re_xc{3} = '0';
    re_yc{3} = '-th2/2';
 re_shape{3} = 'rectangle';
 re_width{3} = 'W';
re_height{3} = 'th2';
re_diel(3) = false;
re_diel_r{3} = '';

  re_name{4} = 'S4';
    re_xc{4} = '0';
    re_yc{4} = 'd/2';
 re_shape{4} = 'rectangle';
 re_width{4} = 'W';
re_height{4} = 'd';
re_diel(4) = true;
re_diel_r{4} = '2.1';

  re_name{5} = 'S5';
    re_xc{5} = '0';
    re_yc{5} = 'd/1.5';
 re_shape{5} = 'rectangle';
 re_width{5} = '1.2*W';
re_height{5} = '2.5*d';
re_diel(5) = true;
re_diel_r{5} = '1';

  re_name{6} = 'S6';
    re_xc{6} = '';
    re_yc{6} = '';
 re_shape{6} = 'rectangle';
 re_width{6} = '';
re_height{6} = '';
re_diel(6) = false;
re_diel_r{6} = '';

  re_name{7} = 'S7';
    re_xc{7} = '';
    re_yc{7} = '';
 re_shape{7} = 'rectangle';
 re_width{7} = '';
re_height{7} = '';
re_diel(7) = false;
re_diel_r{7} = '';

  re_name{8} = 'S8';
    re_xc{8} = '';
    re_yc{8} = '';
 re_shape{8} = 'rectangle';
 re_width{8} = '';
re_height{8} = '';
re_diel(8) = false;
re_diel_r{8} = '';

  re_name{9} = 'S9';
    re_xc{9} = '';
    re_yc{9} = '';
 re_shape{9} = 'rectangle';
 re_width{9} = '';
re_height{9} = '';
re_diel(9) = false;
re_diel_r{9} = '';

  re_name{10} = 'S10';
    re_xc{10} = '';
    re_yc{10} = '';
 re_shape{10} = 'rectangle';
 re_width{10} = '';
re_height{10} = '';
re_diel(10) = false;
re_diel_r{10} = '';


set_formula = 'S1 + S2 + S3 + S3 + S4 + S5';


include_poly = false;

poly_x{1} = 'w';
poly_y{1} = '1';
poly_x{2} = 'd';
poly_y{2} = '1';
poly_x{3} = 'W';
poly_y{3} = '10';
poly_x{4} = 'th1';
poly_y{4} = '0.025';
poly_x{5} = 'th2';
poly_y{5} = '0.2';
poly_x{6} = 'sep';
poly_y{6} = '1';
poly_x{7} = '';
poly_y{7} = '';
poly_x{8} = '';
poly_y{8} = '';


triangle_size = '0.25';
