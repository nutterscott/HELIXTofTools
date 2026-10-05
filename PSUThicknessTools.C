// SC  Following code is a bicubic interpolation around grid points where ToF scintillator thickness measurements were made, with ChatGPT help
// updated, consolidated 20260819 SN 

//In this file:
//void select_four(const double *grid, int n, double q, int idx[4])
//double cubic_lagrange(const double *grid, const double values[],const int idx[4], double q)
//double thickness(int layer, const double *zt, double x, double y)


void select_four(const double *grid, int n, double q, int idx[4])
{
    int first;

    if (q <= grid[1]) {
        first = 0;
    } else if (q >= grid[n - 2]) {
        first = n - 4;
    } else {
        first = 0;
        for (int i = 1; i < n - 2; ++i) {
            if (q < grid[i + 1]) {
                first = i - 1;
                break;
            }
        }
    }

    for (int k = 0; k < 4; ++k)
        idx[k] = first + k;
}

double cubic_lagrange(const double *grid, const double values[],
                      const int idx[4], double q)
{
    double result = 0.0;

    for (int a = 0; a < 4; ++a) {
        double term = values[a];
        for (int b = 0; b < 4; ++b) {
            if (a != b) {
                term *= (q - grid[idx[b]]) /
                        (grid[idx[a]] - grid[idx[b]]);
            }
        }
        result += term;
    }

    return result;
}

double thickness(int layer, const double *zt, double x, double y)
// Accepts paddle or bore 2D layout
{
  if (layer == 0 || layer == 2) //top or bottom
  {
    if (x < 0.0 || x > 20.0 || y < 0.0 || y > 160.0)
        return -1.0;
  }
  else if (layer == 1) //bore
  {
    // SC     the line below was wrong for the bore paddle; the dimensions were for a long paddle
    //    if (x < 0.0 || x > 20.0 || y < 0.0 || y > 160.0)
    if (x < -30.3 || x > 30.3 || y < -30.3 || y > 30.3)
        return -1.0;
  }
  else
  {
    cout << "thickness: Unknown layer: " << layer << endl;
    return -1;
  }

  //top/bottom:
  const int NXtb = 7;
  const int NYtb = 17;
  const double xgridtb[NXtb] = {1., 4., 7., 10., 13., 16., 19.};   // these are fixed for all paddles
  const double ygridtb[NYtb] = {
    1., 10., 20., 30., 40., 50., 60., 70., 80.,
    90., 100., 110., 120., 130., 140., 150., 159.
  };
  //bore:
  const int NXbp = 11;
  const int NYbp = 11;
  const double xgridbp[NXbp] = {-30., -24, -18., -12., -6., 0., 6., 12., 18., 24., 30.};
  const double ygridbp[NYbp] = {-30., -24, -18., -12., -6., 0., 6., 12., 18., 24., 30.};

  int ix[4], iy[4];
  
  int nx = (layer == 1)? NXbp : NXtb;
  int ny = (layer == 1)? NYbp : NYtb;
  const double *xgrid = (layer == 1)? xgridbp : xgridtb;
  const double *ygrid = (layer == 1)? ygridbp : ygridtb;

  //find grid points around (x,y):
  select_four(xgrid, nx, x, ix);
  select_four(ygrid, ny, y, iy);

  double values_y[4];

  for (int j = 0; j < 4; ++j) {
      double values_x[4];
        for (int i = 0; i < 4; ++i)  values_x[i] = zt[ ix[i]*ny + iy[j] ]; // Flattened 2D index formula: row * width + col
      values_y[j] = cubic_lagrange(xgrid, values_x, ix, x);
  }
  

  return cubic_lagrange(ygrid, values_y, iy, y);
}


