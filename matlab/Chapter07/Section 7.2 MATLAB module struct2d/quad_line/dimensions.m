%DIMENSIONS.M - Script used by STRUCT2D.M
  re_name{1} = 'S1';
    re_xc{1} = '0';
    re_yc{1} = 'w+d+th/2';
 re_shape{1} = 'rectangle';
 re_width{1} = 'w';
re_height{1} = 'th';
re_diel(1) = false;
re_diel_r{1} = '';

  re_name{2} = 'S2';
    re_xc{2} = '-(w/2+d+th/2)';
    re_yc{2} = 'w/2';
 re_shape{2} = 'rectangle';
 re_width{2} = 'th';
re_height{2} = 'w';
re_diel(2) = false;
re_diel_r{2} = '';

  re_name{3} = 'S3';
    re_xc{3} = '0';
    re_yc{3} = '-(d+th/2)';
 re_shape{3} = 'rectangle';
 re_width{3} = 'w';
re_height{3} = 'th';
re_diel(3) = false;
re_diel_r{3} = '';

  re_name{4} = 'S4';
    re_xc{4} = 'w/2+d+th/2';
    re_yc{4} = 'w/2';
 re_shape{4} = 'rectangle';
 re_width{4} = 'th';
re_height{4} = 'w';
re_diel(4) = false;
re_diel_r{4} = '';

  re_name{5} = 'S5';
    re_xc{5} = '0';
    re_yc{5} = 'w/2';
 re_shape{5} = 'rectangle';
 re_width{5} = 'w';
re_height{5} = 'w';
re_diel(5) = false;
re_diel_r{5} = '';

  re_name{6} = 'S6';
    re_xc{6} = '0';
    re_yc{6} = 'w+d/2';
 re_shape{6} = 'rectangle';
 re_width{6} = 'w';
re_height{6} = 'd';
re_diel(6) = true;
re_diel_r{6} = '6';

  re_name{7} = 'S7';
    re_xc{7} = '-(w/2+d/2)';
    re_yc{7} = 'w/2';
 re_shape{7} = 'rectangle';
 re_width{7} = 'd';
re_height{7} = 'w';
re_diel(7) = true;
re_diel_r{7} = '6';

  re_name{8} = 'S8';
    re_xc{8} = '0';
    re_yc{8} = '-d/2';
 re_shape{8} = 'rectangle';
 re_width{8} = 'w';
re_height{8} = 'd';
re_diel(8) = true;
re_diel_r{8} = '6';

  re_name{9} = 'S9';
    re_xc{9} = 'w/2+d/2';
    re_yc{9} = 'w/2';
 re_shape{9} = 'rectangle';
 re_width{9} = 'd';
re_height{9} = 'w';
re_diel(9) = true;
re_diel_r{9} = '6';

  re_name{10} = 'S10';
    re_xc{10} = '0';
    re_yc{10} = 'w/2';
 re_shape{10} = 'rectangle';
 re_width{10} = 'w-2*th';
re_height{10} = 'w-2*th';
re_diel(10) = false;
re_diel_r{10} = '';


set_formula = 'S1 + S2 + S3 + S4 +(S5 - S10) + S6 + S7 +S8 + S9 ';


include_poly = false;

poly_x{1} = 'w';
poly_y{1} = '1.27';
poly_x{2} = 'd';
poly_y{2} = '1.27/3';
poly_x{3} = 'th';
poly_y{3} = '0.05';
poly_x{4} = '';
poly_y{4} = '';
poly_x{5} = '';
poly_y{5} = '';
poly_x{6} = '';
poly_y{6} = '';
poly_x{7} = '';
poly_y{7} = '';
poly_x{8} = '';
poly_y{8} = '';


triangle_size = '0.1';
