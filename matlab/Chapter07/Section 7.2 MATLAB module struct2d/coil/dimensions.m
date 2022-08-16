%DIMENSIONS.M - Script used by STRUCT2D.M
  re_name{1} = 'S1';
    re_xc{1} = '0';
    re_yc{1} = '0';
 re_shape{1} = 'ellipse';
 re_width{1} = 'd_out1';
re_height{1} = 'd_out1';
re_diel(1) = false;
re_diel_r{1} = '';

  re_name{2} = 'S2';
    re_xc{2} = '0';
    re_yc{2} = '0';
 re_shape{2} = 'ellipse';
 re_width{2} = 'd_out2';
re_height{2} = 'd_out2';
re_diel(2) = false;
re_diel_r{2} = '';

  re_name{3} = 'S3';
    re_xc{3} = '0';
    re_yc{3} = '0';
 re_shape{3} = 'ellipse';
 re_width{3} = 'd_in1';
re_height{3} = 'd_in1';
re_diel(3) = false;
re_diel_r{3} = '';

  re_name{4} = 'S4';
    re_xc{4} = '0';
    re_yc{4} = '0';
 re_shape{4} = 'ellipse';
 re_width{4} = 'd_in2';
re_height{4} = 'd_in2';
re_diel(4) = false;
re_diel_r{4} = '';

  re_name{5} = 'S5';
    re_xc{5} = '-d_in1/2+cut/2';
    re_yc{5} = '0';
 re_shape{5} = 'rectangle';
 re_width{5} = 'cut';
re_height{5} = 'cut';
re_diel(5) = false;
re_diel_r{5} = '';

  re_name{6} = 'S6';
    re_xc{6} = 'd_in1/2-cut/2';
    re_yc{6} = '0';
 re_shape{6} = 'rectangle';
 re_width{6} = 'cut';
re_height{6} = 'cut';
re_diel(6) = true;
re_diel_r{6} = '40';

  re_name{7} = 'S7';
    re_xc{7} = '0';
    re_yc{7} = '-d_in1/2+cut/2';
 re_shape{7} = 'rectangle';
 re_width{7} = 'cut';
re_height{7} = 'cut';
re_diel(7) = true;
re_diel_r{7} = '1';

  re_name{8} = 'S8';
    re_xc{8} = '0';
    re_yc{8} = 'd_in1/2-cut/2';
 re_shape{8} = 'rectangle';
 re_width{8} = 'cut';
re_height{8} = 'cut';
re_diel(8) = false;
re_diel_r{8} = '';

  re_name{9} = 'S9';
    re_xc{9} = '0';
    re_yc{9} = '0';
 re_shape{9} = 'rectangle';
 re_width{9} = '1.5*d_out1';
re_height{9} = '1.5*d_out1';
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


set_formula = 'S1 - S2 + (S3 - S4) - S5 - S6 - S7 - S8 + (S2 - S3 - S5 - S6 - S7 - S8) + (S4 - S5 - S6 - S7 - S8)';


include_poly = false;

poly_x{1} = 'd_out1';
poly_y{1} = '40';
poly_x{2} = 'd_out2';
poly_y{2} = '39';
poly_x{3} = 'd_in1';
poly_y{3} = '20';
poly_x{4} = 'd_in2';
poly_y{4} = '19';
poly_x{5} = 'cut';
poly_y{5} = '0.7';
poly_x{6} = '';
poly_y{6} = '';
poly_x{7} = '';
poly_y{7} = '';
poly_x{8} = '';
poly_y{8} = '';


triangle_size = '2.0';
