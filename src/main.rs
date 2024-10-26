use std::error::Error;
use std::fs::File;
use std::io::Write;
const DIM: usize = 30;
fn main() -> Result<(), Box<dyn Error>>{
    let nx = 20;
    let ny = 22;
    let mut u = vec![vec![0.0; DIM]; DIM];
    let mut v = vec![vec![0.0; DIM]; DIM];
    let mut un = vec![vec![0.0; DIM]; DIM];
    let mut vn = vec![vec![0.0; DIM]; DIM];
    let mut rho = vec![vec![1.0; DIM]; DIM];
    let mut f = vec![vec![vec![0.0;DIM];DIM]; 9];
    let mut f0 = vec![vec![vec![0.0;DIM];DIM]; 9];
    let mut ftmp = vec![vec![vec![0.0;DIM];DIM]; 9];
    let mut cx = vec![0.0;9];
    let mut cy = vec![0.0;9];
    let gx = 0.00001;
    let gy = 0.0;
    let mut u2;
    let mut norm= 0.0;
    let mut time = 0.0;
    let mut tmp;
    let tau = 0.06;
    let nu = tau/3.0;
    let dx = 1.0;
    let dy = 1.0;
    let dt = 0.1;
    let mut a = vec![vec!['0';DIM];DIM];
  

    // initial condition
    //nu = tau/3.0;

    cx[0] = 0.0; cy[0] = 0.0;
    cx[1] = 1.0; cy[1] = 0.0;  cx[2] = 0.0; cy[2] = 1.0;
    cx[3] =-1.0; cy[3] = 0.0;  cx[4] = 0.0; cy[4] =-1.0;
    cx[5] = 1.0; cy[5] = 1.0;  cx[6] =-1.0; cy[6] = 1.0;
    cx[7] =-1.0; cy[7] =-1.0;  cx[8] = 1.0; cy[8] =-1.0;

    for i in 0..=nx{
        for j in 0..=ny{
            u2 = u[i][j]*u[i][j] + v[i][j]*v[i][j];      
            f0[0][i][j] = rho[i][j]*(1.0 -3.0/2.0*u2)*4.0/9.0;
            for k in 1..=4{
                tmp = cx[k]*u[i][j] + cy[k]*v[i][j];      
                f0[k][i][j] = rho[i][j]*(1.0 +3.0*tmp +9.0/2.0*tmp*tmp -3.0/2.0*u2)/9.0;
            }
            for k in 5..=8{
                tmp = cx[k]*u[i][j] + cy[k]*v[i][j];      
                f0[k][i][j] = rho[i][j]*(1.0 +3.0*tmp +9.0/2.0*tmp*tmp -3.0/2.0*u2)/36.0;
            }
    } }

    for i in 0..=nx{
        for j in 0..=ny{
            for k in 0..=8{
                f[k][i][j] = f0[k][i][j];
    } } }

    for _ in 0..200{
        for _ in 0..1000{
            time = time + dt;

            for i in 0..=nx{
                for j in 0..=ny{
                    un[i][j] = u[i][j]; vn[i][j] = v[i][j];
            } }

            for i in 0..=nx{
                for j in 0..=ny{
                    for k in 0..=8{
                        ftmp[k][i][j] = 0.0;
            } } }

            // collision (SRT)
            for i in 0..=nx{
                for j in 0..=ny{
                u2 = u[i][j]*u[i][j] + v[i][j]*v[i][j];
                f0[0][i][j] = rho[i][j]*(1.0 -3.0/2.0*u2)*4.0/9.0;
                for k in 1..=4{
                    tmp = cx[k]*u[i][j] + cy[k]*v[i][j];      
                    f0[k][i][j] = rho[i][j]*(1.0 +3.0*tmp +9.0/2.0*tmp*tmp -3.0/2.0*u2)/9.0;
                }  
                for k  in 5..=8{
                    tmp = cx[k]*u[i][j] + cy[k]*v[i][j];      
                    f0[k][i][j] = rho[i][j]*(1.0 +3.0*tmp +9.0/2.0*tmp*tmp -3.0/2.0*u2)/36.0;
                }
            } }

            for i in 0..=nx{
                for j in 0..=ny{
                    for k in 0..=8{
                        ftmp[k][i][j] = ftmp[k][i][j] - (f[k][i][j] - f0[k][i][j])/tau;
            } } }

            // force
            for i in 0..=nx{
                for j in 0..=ny{
                    for k in 1..=4{
                        ftmp[k][i][j] = ftmp[k][i][j] + rho[i][j]*(cx[k]*gx + cy[k]*gy)/3.0;
                    }
                    for k in 5..=8{
                        ftmp[k][i][j] = ftmp[k][i][j] + rho[i][j]*(cx[k]*gx + cy[k]*gy)/12.0;
                    }
            } }

            // advection term (FTCS)
            for i in 0..=nx{
                for j in 1..=ny{
                    for k in 0..=8{
                        let mut ip = i as i32+ 1; 
                        if i == nx {ip =  0;}
                        let mut im = i as i32 - 1; 
                        if i ==  0 {im = nx as i32;}
                        ftmp[k][i][j] = ftmp[k][i][j] - cx[k]*(f[k][ip as usize][j  ] - f[k][im as usize][j  ])/dx*0.5 
                                    - cy[k]*(f[k][i ][j+1] - f[k][i ][j-1])/dy*0.5; 
            } } }

            // boundary condition (mesoscale, extrapolation)
            for i in 0..=nx{
                for k in 0..=8{
                    ftmp[k][i][ 0] = 2.0*ftmp[k][i][   1] - ftmp[k][i][   2];
                    ftmp[k][i][ny] = 2.0*ftmp[k][i][ny-1] - ftmp[k][i][ny-2];
            } }

            // propagation (FTCS)
            for i in 0..=nx{
                for j in 0..=ny{
                    for k in 0..=8{
                        f[k][i][j] = f[k][i][j] + ftmp[k][i][j]*dt;
            } } }

            // physics
            for i in 0..=nx{
                for j in 0..=ny{
                    rho[i][j] = f[0][i][j]; u[i][j] = 0.0; v[i][j] = 0.0;
                for k in 1..=8{
                    rho[i][j] = rho[i][j] + f[k][i][j];
                    u[i][j] =   u[i][j] + f[k][i][j]*cx[k];
                    v[i][j] =   v[i][j] + f[k][i][j]*cy[k];
                } 
                u[i][j] = u[i][j]/rho[i][j];
                v[i][j] = v[i][j]/rho[i][j];
            } }

            // boundary condition (macroscale)
            for i in 0..=nx{
                u[i][   1] = 0.0; v[i][   1] = 0.0;
                u[i][ny-1] = 0.0; v[i][ny-1] = 0.0;
            }

            norm = 0.0;
            for i in 0..=nx{
                for j in 0..=ny{
                    tmp = ( f64::powf(u[i][j]-un[i][j],2.0) + f64::powf(v[i][j]-vn[i][j], 2.0)).sqrt();
                    if tmp > norm{ 
                        norm = tmp;
                    }
            } }

  } //loop2
  println!("Time = {:10.2}, Norm = {:10.2}", time, norm);
  println!("Umax = {:8.6e}({:6.4e})", u[nx/2][ny/2], gx/8.0/nu*(ny as f64 - 2.0)*(ny as f64 - 2.0));
  
  for i in 0..=nx{
    for j in 0..=ny{
        a[i][j]='0';
  } }

  for j in 0..=ny{
    let nb = u[nx/2][j]/(gx/8.0/nu*(ny as f64-2.0)*(ny as f64 - 2.0))*20.0;
    let mut i = 0;
    while i as f64 <= nb{
        a[i][j] = '-';
        i+=1;
    }
  }

  for j in 0..=ny{
    for i in 0..=nx{
        print!("{}", a[i][j]);
    }
    println!("");
  }

  if norm < 0.0000000001 && time > 10000.0{ 
    let mut file = File::create("data")?;
    for j in 0..=ny{
        writeln!(file, "{:10.8e}", u[nx/2][j]/(gx/8.0/nu*(ny as f64-2.0)*(ny as f64 - 2.0)))?;
    }
    
  }

  } //loop1

  Ok(())

}