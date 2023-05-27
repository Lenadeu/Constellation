using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO; 

namespace ConstellationPrj
{
    public partial class MainForm : Form
    {
        private static ToDegOrRad _Conv = new ToDegOrRad();

        public MainForm()
        {
            InitializeComponent();
        }

        private void Go_Click(object sender, EventArgs e)
        {
            ObjectState Init;
            // Initials 
            Double Alfa0, Betta0, Omega_x0, Omega_y0, Omega_z0, Theta0, Psi0, Gamma0, X0, Y0, Z0, m0;
            Double T, dt, t0;
            Boolean SchemaType = false;
            Boolean ControlType = false;
            //--------------------------------------------------------------
            Alfa0 = _Conv.ToRad(1);
            Betta0 = 0;
            Omega_x0 = _Conv.ToRad(0);
            Omega_y0 = _Conv.ToRad(0);
            Omega_z0 = _Conv.ToRad(0);
            Theta0 = _Conv.ToRad(4);
            Psi0 = _Conv.ToRad(2);
            Gamma0 = _Conv.ToRad(0);
            X0 = 0;
            Y0 = 2458.2;
            Z0 = 0;
            m0 = 2161;
            Init = new ObjectState(Alfa0, Betta0, Omega_x0, Omega_y0, Omega_z0, Theta0, Psi0, Gamma0, X0, Y0, Z0, m0);
            //-------------------------------------------------------------
            // Total time of computation 
            T = 2;
            // Шаг расчета 
            dt = 0.0001;
            // Starting value of time  
            t0 = 11;
            // Aircraft transient computation  
            Transient Process = new Transient(T, dt, t0, Init, SchemaType, ControlType, trjCompProgressBar);
            
            Process.CalculateTransient();

        }

      

    
    }
    // Transient 
    class Transient
    {
        private Double _T; // Total time of computation 
        private Double _dt; // Integration step 
        private Double _t0; // Starting time  
        private Int32 _N; // Number of iterations        

        private Object _Obj; // Aircraft as an object of control  
        private static ToDegOrRad _Conv = new ToDegOrRad();
        ProgressBar _IntegrationProcessBar;

        
        public Transient(Double T, Double dt, Double t0, ObjectState Initial, Boolean SchemaType, Boolean ControlType, ProgressBar IntegrationProcessBar)
        {
            _T = T;
            _dt = dt;
            _t0 = t0;
            _N = Convert.ToInt32(Math.Floor(_T / _dt)); 
            _Obj = new Object(Initial, SchemaType);
            _IntegrationProcessBar = IntegrationProcessBar;
            
        }

        //  
        public void CalculateTransient()
        {
            Double t = 0; // Current time 
            Int32 k = 0; // Number of iteration   
            StreamWriter Output = new StreamWriter("Output" + "_" + DateTime.Now.Date.Day.ToString() + "_" + DateTime.Now.Hour.ToString() + "_" + DateTime.Now.Minute.ToString() + ".txt");

            // Time stamp for writing data into output file 
            Double tp = 0;
            // Step of writing data into output file 
            Double dtp = 0.0005;
    

            // Initialization of progress bar 
            _IntegrationProcessBar.Minimum = 0;
            _IntegrationProcessBar.Maximum = _N;
            _IntegrationProcessBar.Step = 1;        
          

            while (k < _N)
            {

                // One step forward for progress bar 
                _IntegrationProcessBar.PerformStep();


                t = _t0 + k * _dt;

                try
                {
                   
                    // One step of intgration of aircraft equations                  
                    _Obj.GetNextState(_dt, t);
                   
                    // Output into files 
                    if (t >= tp)
                    {
                        tp += dtp;
                        Output.Write(String.Format("{0,-15:F4}", t));
                       
                        Output.Write(String.Format("{0,-15:F4}", _Obj.CurrentState.X));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.CurrentState.Y));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.CurrentState.Z));
                        Output.Write(String.Format("{0,-15:F4}", _Conv.ToDeg(_Obj.CurrentState.Theta)));
                        Output.Write(String.Format("{0,-15:F4}", _Conv.ToDeg(_Obj.CurrentState.Psi)));
                        Output.Write(String.Format("{0,-15:F10}", _Conv.ToDeg(_Obj.CurrentState.Gamma)));
                        Output.Write(String.Format("{0,-15:F6}", _Obj.Aero.GammaCalc));

                        Output.Write(String.Format("{0,-15:F4}", _Conv.ToDeg(_Obj.CurrentState.Omega_x)));
                        Output.Write(String.Format("{0,-15:F4}", _Conv.ToDeg(_Obj.CurrentState.Omega_y)));
                        Output.Write(String.Format("{0,-15:F4}", _Conv.ToDeg(_Obj.CurrentState.Omega_z)));

                        Output.Write(String.Format("{0,-15:F4}", _Obj.Cy));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.Cz));

                        Output.Write(String.Format("{0,-15:F4}", _Obj.V));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.q));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.ro));

                        Output.Write(String.Format("{0,-15:F4}", _Obj.Mx));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.My));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.Mz));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.Aero.mx));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.Aero.my));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.Aero.mz));

                        Output.Write(String.Format("{0,-15:F4}", _Obj.nx));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.ny));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.nz));

                        Output.Write(String.Format("{0,-15:F4}", _Obj.Jx));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.Jy));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.Jz));

                        Output.Write(String.Format("{0,-15:F4}", _Obj.dXc));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.CentrF.Xc));

                      


                        Output.Write(String.Format("{0,-15:F4}", _Obj.M));

                        Output.Write(String.Format("{0,-15:F4}", _Obj.Aero.Cy0));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.Aero.dCyd1));

                        Output.Write(String.Format("{0,-15:F4}", _Obj.Aero.dCyd2));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.Aero.dCyd3));
                        Output.Write(String.Format("{0,-15:F4}", _Obj.Aero.dCyd4));
                        Output.Write(String.Format("{0,-15:F6}", _Obj.Aero.Alfa_s));
                    

                        Output.WriteLine();

                    }
                }
                catch
                {

                    Output.Close();
                    MessageBox.Show("An error in the integration loop!!!");
                    break;
                }

                k++;

            }
            MessageBox.Show("Computation is fnished! Results are written into output text file");
            Output.Close();

        }


    }

    // Degress to radians transformation (and back)
    class ToDegOrRad
    {

        public Double ToRad(Double value)
        {
            return value * Math.PI / 180;

        }
        public Double ToDeg(Double value)
        {
            return value * 180 / Math.PI;
        }
    }
    // Interpolation 
    class Interpolation
    {
        public Double Interp1D(Double[] x, Double[] y, Double xi)
        {
            Int16 i;
            Int16 n = Convert.ToInt16(x.Length);
            Double h, e;
            //  Find first index of element of array x which is more than xi
            if (xi < x[0]) //  we are out of left boundary 
            {
                i = 0; // 
            }
            else
            {
                if (xi > x[n - 1]) //  we are out of right boundary
                {
                    i = Convert.ToInt16(n - 2); // 
                }
                else
                {
                    for (i = 0; i < n - 2; i++)
                    {
                        if (xi >= x[i] && xi <= x[i + 1]) // xi in between of x[i] and x[i+1] 
                        {
                            break;
                        }
                    }
                }
            }
            // Step regarding to x 
            h = x[i + 1] - x[i];

            e = (xi - x[i]) / h; // Relation of the distance between xi and x[i] and the distance between x[i] and x[i+1]
            return y[i] * (1 - e) + y[i + 1] * e;

        }

        public Double Interp2D(Double[] x, Double[] y, Double[,] z, Double xi, Double yi)
        {

            Int16 i, j, ny, nx;
            Double h, e, fx1, fx2;
            nx = Convert.ToInt16(x.Length);
            ny = Convert.ToInt16(y.Length);
            Double[] x_sr = new Double[ny];
            if (xi < x[0])
            {
                i = 0;
            }
            else
            {
                if (xi > x[nx - 1])
                {
                    i = Convert.ToInt16(nx - 2);
                }
                else
                {
                    for (i = 0; i < nx - 2; i++)
                    {
                        if ((xi >= x[i]) && (xi <= x[i + 1]))
                        {
                            break;
                        }
                    }
                }
            }
            for (j = 0; j < ny; j++)
            {
                x_sr[j] = z[i, j];
            }
            fx1 = Interp1D(y, x_sr, yi);
            for (j = 0; j < ny; j++)
            {
                x_sr[j] = z[i + 1, j];
            }
            fx2 = Interp1D(y, x_sr, yi);
            h = x[i + 1] - x[i];
            e = (xi - x[i]) / h;
            return fx1 * (1 - e) + fx2 * e;

        }
        public Double Interp3D(Double[] x, Double[] y, Double[] z, Double[, ,] u, Double xi, Double yi, Double zi)
        {
            Int16 i, j, k, nx, ny, nz;
            Double fxy1, fxy2, h, e;
            nx = Convert.ToInt16(x.Length);
            ny = Convert.ToInt16(y.Length);
            nz = Convert.ToInt16(z.Length);
            Double[,] xy_sr = new Double[ny, nz];
            if (xi < x[0])
            {
                i = 0;
            }
            else
            {
                if (xi > x[nx - 1])
                {
                    i = Convert.ToInt16(nx - 2);
                }
                else
                {
                    for (i = 0; i < nx - 2; i++)
                    {
                        if (xi >= x[i] && xi <= x[i + 1])
                        {
                            break;
                        }
                    }
                }
            }
            for (j = 0; j < ny; j++)
            {
                for (k = 0; k < nz; k++)
                {
                    xy_sr[j, k] = u[i, j, k];

                }
            }
            fxy1 = Interp2D(y, z, xy_sr, yi, zi);
            for (j = 0; j < ny; j++)
            {
                for (k = 0; k < nz; k++)
                {
                    xy_sr[j, k] = u[i + 1, j, k];
                }
            }
            fxy2 = Interp2D(y, z, xy_sr, yi, zi);
            h = x[i + 1] - x[i];
            e = (xi - x[i]) / h;
            return fxy1 * (1 - e) + fxy2 * e;

        }

        public Double Interp4D(Double[] x, Double[] y, Double[] z, Double[] u, Double[, , ,] p, Double xi, Double yi, Double zi, Double ui)
        {
            Int16 i, j, k, t, nx, ny, nz, nu;
            Double fxyz1, fxyz2, h, e;
            nx = Convert.ToInt16(x.Length);
            ny = Convert.ToInt16(y.Length);
            nz = Convert.ToInt16(z.Length);
            nu = Convert.ToInt16(u.Length);
            Double[, ,] xyz_sr = new Double[ny, nz, nu];
            if (xi < x[0])
            {
                i = 0;
            }
            else
            {
                if (xi > x[nx - 1])
                {
                    i = Convert.ToInt16(nx - 2);
                }
                else
                {
                    for (i = 0; i < nx - 2; i++)
                    {
                        if (xi >= x[i] && xi <= x[i + 1])
                        {
                            break;
                        }
                    }
                }
            }
            for (t = 0; t < ny; t++)
            {
                for (j = 0; j < nz; j++)
                {
                    for (k = 0; k < nu; k++)
                    {
                        xyz_sr[t, j, k] = p[i, t, j, k];
                    }
                }
            }

            fxyz1 = Interp3D(y, z, u, xyz_sr, yi, zi, ui);
            for (t = 0; t < ny; t++)
            {
                for (j = 0; j < nz; j++)
                {
                    for (k = 0; k < nu; k++)
                    {
                        xyz_sr[t, j, k] = p[i + 1, t, j, k];
                    }
                }
            }

            fxyz2 = Interp3D(y, z, u, xyz_sr, yi, zi, ui);
            h = x[i + 1] - x[i];
            e = (xi - x[i]) / h;
            return fxyz1 * (1 - e) + fxyz2 * e;


        }




    }
    // State of object of control 
    class ObjectState
    { 
        private Double _V; // Скорость  Velicity 
        private Double _Theta_b; // Vertical path angle  
        private Double _Psi_b; //  Horizontal path angle 
        private Double _Omega_x; //  Angular velocity regarding to x  
        private Double _Omega_y; //  Angular velocity regarding to y
        private Double _Omega_z; // Angular velocity regarding to z
        private Double _Theta; // Pitch angle 
        private Double _Psi; // Yaw angle 
        private Double _Gamma; // Roll angle 
        private Double _X; // Horizontal coordinate 
        private Double _Y; // Vertical coordinate 
        private Double _Z; // Lateral coordinate 
        
        public Double V
        {
            get
            {
                return _V;
            }
            set
            {
                _V = value;
            }
        }
        public Double Theta_b
        {
            get
            {
                return _Theta_b;
            }
            set
            {
                _Theta_b = value;
            }
        }
        public Double Psi_b
        {
            get
            {
                return _Psi_b;
            }
            set
            {
                _Psi_b = value;
            }
        }
        public Double Omega_x
        {
            get
            {
                return _Omega_x;
            }
            set
            {
                _Omega_x = value;
            }
        }
        public Double Omega_y
        {
            get
            {
                return _Omega_y;
            }
            set
            {
                _Omega_y = value;
            }
        }
        public Double Omega_z
        {
            get
            {
                return _Omega_z;
            }
            set
            {
                _Omega_z = value;
            }
        }
        public Double Theta
        {
            get
            {
                return _Theta;
            }
            set
            {
                _Theta = value;
            }
        }
        public Double Psi
        {
            get
            {
                return _Psi;
            }
            set
            {
                _Psi = value;
            }
        }
        public Double Gamma
        {
            get
            {
                return _Gamma;
            }
            set
            {
                _Gamma = value;
            }
        }
        public Double X
        {
            get
            {
                return _X;
            }
            set
            {
                _X = value;
            }
        }
        public Double Y
        {
            get
            {
                return _Y;
            }
            set
            {
                _Y = value;
            }
        }
        public Double Z
        {
            get
            {
                return _Z;
            }
            set
            {
                _Z = value;
            }
        }
       
        public ObjectState()
        {
            _V = 0;
            _Theta_b = 0;
            _Psi_b = 0;
            _Omega_x = 0;
            _Omega_y = 0;
            _Omega_z = 0;
            _Theta = 0;
            _Psi = 0;
            _Gamma = 0;
            _X = 0;
            _Y = 0;
            _Z = 0;
        }

        public ObjectState(Double V, Double Theta_b, Double Psi_b, Double Omega_x, Double Omega_y, Double Omega_z, Double Theta, Double Psi, Double Gamma, Double X, Double Y, Double Z)
        {
            _V = V;
            _Theta_b = Theta_b;
            _Psi_b = Psi_b;
            _Omega_x = Omega_x;
            _Omega_y = Omega_y;
            _Omega_z = Omega_z;
            _Theta = Theta;
            _Psi = Psi;
            _Gamma = Gamma;
            _X = X;
            _Y = Y;
            _Z = Z;
        }
        public ObjectState(ObjectState InitObject)
        {

            _V = InitObject.V;
            _Theta_b = InitObject.Theta_b;
            _Psi_b = InitObject.Psi_b;
            _Omega_x = InitObject.Omega_x;
            _Omega_y = InitObject.Omega_y;
            _Omega_z = InitObject.Omega_z;
            _Theta = InitObject.Theta;
            _Psi = InitObject.Psi;
            _Gamma = InitObject._Gamma;
            _X = InitObject.X;
            _Y = InitObject.Y;
            _Z = InitObject.Z;
            
        }
        public void SetObjectState(ObjectState Object)
        {
            _V = Object.V;
            _Theta_b = Object.Theta_b;
            _Psi_b = Object.Psi_b;
            _Omega_x = Object.Omega_x;
            _Omega_y = Object.Omega_y;
            _Omega_z = Object.Omega_z;
            _Theta = Object.Theta;
            _Psi = Object.Psi;
            _Gamma = Object._Gamma;
            _X = Object.X;
            _Y = Object.Y;
            _Z = Object.Z;
            
        }

    }

    // Object 
    class Object
    {
        private ObjectState _CurSt; // Current state of control object  
        private ObjectState _PrevSt; // Previous state of control object 
       
        private Atmosphere _Atm; // Atmosphere 
       
        private Aerodynamics _Aero; // Aerodynamics 
        private InertialFeatures _InertF; // Inertia moments  
        private CentringFeatures _CentrF; // Centring characteristics
        private MidshipCharacters _MidCh; // Geometric characteristics of the midsection 
        // 
        private Double _Cx,_Cy, _Cz, _Jx, _Jy, _Jz, _V, _M, _Mx, _My, _Mz, _q, _ro, _dXc, _S, _D,_m,_alfa,_betta,_gamma_c;
        // Trajectory roll angle 
        public Double _Gammag;
        // g 
        private const Double _g = 9.81;
        // Overloads 
        private Double _nx, _ny, _nz;
        // 
        private static ToDegOrRad _Conv = new ToDegOrRad();
              
        public Double Mx
        {
            get
            {
                return _Mx;
            }
        }
        public Double My
        {
            get
            {
                return _My;
            }
        }
        public Double Mz
        {
            get
            {
                return _Mz;
            }
        }
        public Double q
        {
            get
            {
                return _q;
            }
        }

        public Double ro
        {
            get
            {
                return _ro;
            }
        }

        public Double V
        {
            get
            {
                return _V;
            }
        }

       
        public Aerodynamics Aero
        {
            get
            {
                return _Aero;
            }
        }
        public Double Jx
        {
            get
            {
                return _Jx;
            }
        }
        public Double Jy
        {
            get
            {
                return _Jy;
            }
        }
        public Double Jz
        {
            get
            {
                return _Jz;
            }
        }
        public Double dXc
        {
            get
            {
                return _dXc;
            }
        }
        public CentringFeatures CentrF
        {
            get
            {
                return _CentrF;
            }
        }
        public ObjectState CurrentState
        {
            get
            {
                return _CurSt;
            }
        }
        public Double nx
        {
            get
            {
                return _nx;
            }
        }
        public Double ny
        {
            get
            {
                return _ny;
            }
        }
        public Double nz
        {
            get
            {
                return _nz;
            }
        }
        public Double Cx
        {
            get
            {
                return _Cy;
            }
        }

        public Double Cy
        {
            get
            {
                return _Cy;
            }
        }
        public Double Cz
        {
            get
            {
                return _Cz;
            }
        }
        public Double M
        {
            get
            {
                return _M;
            }
        }
     

        // Input arguments: Initials and type of wings schematics   
        public Object(ObjectState Initial, Boolean SchemaType)
        {
            Int32 i;
            _PrevSt = new ObjectState(Initial);
            _CurSt = new ObjectState();
           
            _Atm = new Atmosphere();
           
            _Aero = new Aerodynamics("Cy0.txt", "dCy.txt", "Cz0.txt", "dCz.txt", "mx0.txt", "dmx.txt", "my0.txt", "dmy.txt", "mz0.txt", "dmz.txt", SchemaType);
            _InertF = new InertialFeatures();
            _CentrF = new CentringFeatures();
            _MidCh = new MidshipCharacters();
            _nx = 0;
            _ny = 0;
            _nz = 0;

          
        }
      
        // One step of integration of aircraft equations 
        public void GetNextState(Double dt, Double t)
        {
            Int16 i;
            Double[] DeltasT = new Double[4];

            
            // Radians to grads transition  
            for (i = 0; i < 4; i++)
            {
                DeltasT[i] = _Conv.ToDeg(DeltasT[i]);
            }
            // Mass
            _m = 800;
            // Mach number  
            _M = _Atm.GetM(_PrevSt.Y, _PrevSt.V);
            // Air density  
            _ro = _Atm.Get_ro(_PrevSt.Y);
            // Speed pressure   
            _q = _ro * Math.Pow(_PrevSt.V, 2) / 2;
            //   
            _dXc = _CentrF.Get_dXc(_m);
            // Square of midel 
            _S = _MidCh.S;
            // Dameter of midel  
            _D = _MidCh.D;

            _alfa = _PrevSt.Theta - _PrevSt.Theta_b;
            _betta = Math.Cos(_PrevSt.Theta_b) * _PrevSt.Psi - Math.Cos(_PrevSt.Theta) * _PrevSt.Psi_b + _alfa * _PrevSt.Gamma;

            // Aerodynamic coefficients 
            _Aero.CalculateAerodynamics(_Conv.ToDeg(_alfa), _Conv.ToDeg(_betta), _M, DeltasT, _dXc, _D,_S, _q);
            _Cx = _Aero.Cx;
            _Cy = _Aero.Cy;
            _Cz = _Aero.Cz;
            _Mx = _Aero.Mx;
            _My = _Aero.My;
            _Mz = _Aero.Mz;         
            // Moments of inertia  
            _Jx = _InertF.Get_Jx(_m);
            _Jy = _InertF.Get_Jy(_m);
            _Jz = _InertF.Get_Jz(_m);

            _gamma_c=Math.Tan(_PrevSt.Theta_b)*_betta+Math.Cos(_PrevSt.Theta)/Math.Cos(_PrevSt.Theta_b)*_PrevSt.Gamma;
            
            // Integration 

            _CurSt.V = _PrevSt.V + dt * ( 1 / _m * (-_Cx * _q * _S - _m * _g * Math.Sin(_PrevSt.Theta_b)) );

            _CurSt.Theta_b = _PrevSt.Theta_b + dt * (1 / (_m * _PrevSt.V) * (-_Cz * _q * _S * Math.Sin(_gamma_c) + _Cy * _q * _S * Math.Cos(_gamma_c) - _m * _g * Math.Cos(_PrevSt.Theta_b)));

            _CurSt.Psi_b = _PrevSt.Psi_b + dt * (-1 / (_m * _PrevSt.V) * (_Cy *_q * _S * Math.Sin(_gamma_c) + _Cz* _q * _S * Math.Cos(_gamma_c))) ;

            _CurSt.Omega_x = _PrevSt.Omega_x + dt * (1/_Jx*(_Mx - (_Jz - _Jy) * _PrevSt.Omega_y * _PrevSt.Omega_z));

            _CurSt.Omega_y = _PrevSt.Omega_y + dt * (1 / _Jy * (_My - (_Jx - _Jz) * _PrevSt.Omega_z * _PrevSt.Omega_x));

            _CurSt.Omega_z = _PrevSt.Omega_z + dt * (1 / _Jz * (_Mz - (_Jy - _Jx) * _PrevSt.Omega_x * _PrevSt.Omega_y));

            _CurSt.Psi = _PrevSt.Psi + dt * ((_PrevSt.Omega_y * Math.Cos(_PrevSt.Gamma) - _PrevSt.Omega_z * Math.Sin(_PrevSt.Gamma)) / Math.Cos(_PrevSt.Theta));

            _CurSt.Theta = _PrevSt.Theta + dt * ((_PrevSt.Omega_y * Math.Sin(_PrevSt.Gamma) + _PrevSt.Omega_z * Math.Cos(_PrevSt.Gamma)));

            _CurSt.Gamma = _PrevSt.Gamma + dt * ( _PrevSt.Omega_x - (_PrevSt.Omega_y * Math.Cos(_PrevSt.Gamma) - _PrevSt.Omega_z * Math.Sin(_PrevSt.Gamma)) * Math.Tan(_PrevSt.Theta) );

            _CurSt.X = _PrevSt.X + dt * (_PrevSt.V * Math.Cos(_PrevSt.Theta_b) * Math.Cos(_PrevSt.Psi_b));

            _CurSt.Y = _PrevSt.Y + dt * (_PrevSt.V * Math.Sin(_PrevSt.Theta_b));

            _CurSt.Z = _PrevSt.Z + dt * (-_PrevSt.V * Math.Sin(_PrevSt.Psi_b) * Math.Cos(_PrevSt.Theta_b));

            // Assigning previous state to current state   
            _PrevSt.SetObjectState(_CurSt);
        }

       
    }

    // Geometric characteristics of the midsection
    class MidshipCharacters
    {
        // Midsection Diameter, m
      private const Double _D = 0.665;
      // Midsection square , m^2 
      private const Double _S = 0.347;

      public Double D      
      {
          get
          {
              return _D; 
          }
      }
      public Double S
      {
          get
          {
              return _S;
          }
      }
    }



   // Aerodynamics 
    class Aerodynamics
    {
        
        // Table for calculating the normal force coefficient when the controls are not tilted
        private FileData2D _Cy0Tabl;
        // Table for calculating the incremental normal force coefficient due to deviation of one rudder
        private FileData4D _dCyTabl;
        // Normal force coefficient with controls not tilted
        private Double _CyT;
        // Normal force coefficient for undeflected controls in the associated marshalling stage coordinate system
        private Double _Cy;

        // Table for calculating the lateral force coefficient when the controls are not tilted
        private FileData3D _Cz0Tabl;
        // Table for calculating the incremental lateral force coefficient due to deviation of one rudder
        private FileData4D _dCzTabl;
        // Lateral force coefficient when controls are not tilted
        private Double _CzT;
        // Lateral force coefficient for undeflected controls in the associated marshalling stage coordinate system
        private Double _Cz;

        private Double _Cx;


        // Table for calculating the yaw moment coefficient when the controls are not tilted
        private FileData3D _my0Tabl;
        // Table for calculating the incremental yaw moment coefficient due to deviation of one rudder
        private FileData4D _dmyTabl;
        // Yaw moment coefficient 
        private Double _myT;
        // Yaw moment coefficient in the associated marshalling stage coordinate system
        private Double _my;

        //  Table for calculating the pitching moment coefficient when the controls are not tilted
        private FileData2D _mz0Tabl;
        // Table for calculating the increase in pitch factor due to deviation of one rudder
        private FileData4D _dmzTabl;
        // Pitch moment coefficient
        private Double _mzT;
        // Pitching moment coefficient in the associated marshalling stage coordinate system
        private Double _mz;

        // Table for calculating the roll moment coefficient when the controls are not tilted
        private FileData3D _mx0Tabl;
        // Table for calculating the incremental roll moment coefficient due to deviation of one rudder
        private FileData4D _dmxTabl;
        // Roll moment coefficient 
        private Double _mxT;
        // Rolling moment coefficient in the associated marshalling stage coordinate system
        private Double _mx;

        
        private Interpolation _Interp;
        
        private static ToDegOrRad _Conv = new ToDegOrRad();
      
        private Double _Cy0, _dCyd1, _dCyd2, _dCyd3, _dCyd4;
        private Double _Cz0, _dCzd1, _dCzd2, _dCzd3, _dCzd4;
        private Double _my0, _dmyd1, _dmyd2, _dmyd3, _dmyd4;
        private Double _mz0, _dmzd1, _dmzd2, _dmzd3, _dmzd4;
        private Double _mx0, _dmxd1, _dmxd2, _dmxd3, _dmxd4;
        private Double _Mx, _My, _Mz;
        private Double _Yact, _Zact;
        private Double _CyTact,_CzTact;

        private Double _GammaCalc;
        private Double _Alfa_s;
        private Double[] _GammaActs = new Double[4];
        //  Rudder orientation schematics. 0 - "+" schematics, 1 - "x" schematics 
        private Boolean _SchemaType; 

        
        // Accuracy 
        private const Double eps = 0.0000001;

        public Double Cy0
        {
            get
            {
                return _Cy0; 

            }
        }

        public Double dCyd1
        {
            get
            {
                return _dCyd1;

            }
        }
        public Double dCyd2
        {
            get
            {
                return _dCyd2;

            }
        }
        public Double dCyd3
        {
            get
            {
                return _dCyd3;

            }
        }
        public Double dCyd4
        {
            get
            {
                return _dCyd4;

            }
        }

        public Double GammaCalc
        {
            get
            {
                return _GammaCalc;
 
            }
        }
        public Double Alfa_s
        {
            get
            {
                return _Alfa_s;

            }
        }
        public Double mx
        {
            get
            {
                return _mx;
            }
        }
        public Double my
        {
            get
            {
                return _my;
            }
        }
        public Double mz
        {
            get
            {
                return _mz;
            }
        }
        public Double Mx
        {
            get
            {
                return _Mx;
            }
        }
        public Double My
        {
            get
            {
                return _My;
            }
        }
        public Double Mz
        {
            get
            {
                return _Mz;
            }
        }
        public Double Cy
        {
            get
            {
                return _Cy;
            }
        }
        public Double Cz
        {
            get
            {
                return _Cz;
            }
        }
        public Double Cx
        {
            get
            {
                return _Cz;
            }
        }
        public Double Yact
        {
            get
            {
                return _Yact; 
            }
        }
        public Double Zact
        {
            get
            {
                return _Zact;
            }
        }


        public Aerodynamics(String Cy0_f_name, String dCy_f_name, String Cz0_f_name, String dCz_f_name, String mx0_f_name, String dmx_f_name, String my0_f_name, String dmy_f_name, String mz0_f_name, String dmz_f_name, Boolean SchemaType)
        {
            // Object for reading data from aerodynamic tables 
            ReadData RD;
            RD = new ReadData();

            _Cy0Tabl = RD.ReadDataFrom2DTable(Cy0_f_name);
            _dCyTabl = RD.ReadDataFrom4DTable(dCy_f_name);

            _Cz0Tabl = RD.ReadDataFrom3DTable(Cz0_f_name);
            _dCzTabl = RD.ReadDataFrom4DTable(dCz_f_name);

            _mx0Tabl = RD.ReadDataFrom3DTable(mx0_f_name);
            _dmxTabl = RD.ReadDataFrom4DTable(dmx_f_name);

            _my0Tabl = RD.ReadDataFrom3DTable(my0_f_name);
            _dmyTabl = RD.ReadDataFrom4DTable(dmy_f_name);

            _mz0Tabl = RD.ReadDataFrom2DTable(mz0_f_name);
            _dmzTabl = RD.ReadDataFrom4DTable(dmz_f_name);

            _Interp = new Interpolation();
            _SchemaType = SchemaType;
            
        }
        public void CalculateAerodynamics(Double Alfa, Double Betta,Double M, Double[] Deltas, Double dXc,Double D, Double S, Double q)
        {
        
            Int32 i; 
            _Alfa_s = GetAlfa_s(Alfa, Betta);
            _GammaCalc = GetGammaCalc(Alfa, Betta, _Alfa_s);
            
            if (_GammaCalc<0)
            {
                _GammaCalc += 360;
            }
            for (i = 0; i < 4; i++)
            {
                _GammaActs[i] = GetGammaAct(i+1);
            }
            _CyT = GetCyT(M, Deltas);
            _CzT = GetCzT(M, Deltas);
            // Normal force coefficient in the coordinate system associated with the product
            _Cy = _CyT* Math.Cos(_Conv.ToRad(_GammaCalc)) + _CzT * Math.Sin(_Conv.ToRad(_GammaCalc));
            // Lateral force coefficient in the coordinate system associated with the product 
            _Cz = _CzT*Math.Cos(_Conv.ToRad(_GammaCalc))-_CyT*Math.Sin(_Conv.ToRad(_GammaCalc));
            _myT = GetmyT(M, Deltas);
            _mzT = GetmzT(M, Deltas);
            // Conversion to the associated aircraft coordinate system
            _my = _myT*Math.Cos(_Conv.ToRad(_GammaCalc))+_mzT*Math.Sin(_Conv.ToRad(_GammaCalc));
            // Recalculating to the new alignment 
            _my=_my-(dXc / D) * _Cz;
            // Conversion to the coordinate system associated with the product
            _mz = _mzT * Math.Cos(_Conv.ToRad(_GammaCalc)) -_myT* Math.Sin(_Conv.ToRad(_GammaCalc));
            // Recalculating to the new alignment 
            _mz = _mz + (dXc / D) * _Cy;
            // Conversion to the associated aircraft coordinate system
            _mx = GetmxT(M, Deltas);
            // Rolling moment in the coordinate system associated with the marshalling stage
             _Mx = _mx* q * S * D;
             // Momentum of yaw in the coordinate system associated with the aircraft 
            _My = _my * q * S * D;
            // Pitching torque in the associated marshalling stage coordinate system
            _Mz = _mz * q * S * D;
            _CyTact = _dCyd1 + _dCyd2 + _dCyd3 + _dCyd4;
            _CzTact = _dCzd1 + _dCzd2 + _dCzd3 + _dCzd4;
            _Yact = (_CyTact* Math.Cos(_Conv.ToRad(_GammaCalc))+_CzTact*Math.Sin(_Conv.ToRad(_GammaCalc)))* q * S;
            _Zact = (_CzTact * Math.Cos(_Conv.ToRad(_GammaCalc)) - _CyTact * Math.Sin(_Conv.ToRad(_GammaCalc))) * q * S;
        
        }
        // Calculation of spatial angle of attack through angles of attack and glide. Entry and exit in degrees!!! 
        private Double GetAlfa_s(Double Alfa, Double Betta)
        {
            Double Alfa_s;
            Alfa_s = Math.Acos(Math.Cos(_Conv.ToRad(Alfa)) * Math.Cos(_Conv.ToRad(Betta)));        
            Alfa_s = _Conv.ToDeg(Alfa_s);
            return Alfa_s;
        }

        // Calculation of the angle of rotation around the aircraft axis. Entry and exit in degrees
        private Double GetGammaCalc(Double Alfa, Double Betta, Double Alfa_s)
        {
            Double Gamma, SinGamma, CosGamma;

            if (Math.Abs(Math.Sin(_Conv.ToRad(Alfa_s))) < eps) // alfa=0 and betta=0, gamma = 0 
            {
                SinGamma = 0;
                CosGamma = 1;
            }
            else 
            {
                SinGamma = Math.Sin(_Conv.ToRad(Betta)) / Math.Sin(_Conv.ToRad(Alfa_s));
                CosGamma = Math.Sin(_Conv.ToRad(Alfa)) * Math.Cos(_Conv.ToRad(Betta)) / Math.Sin(_Conv.ToRad(Alfa_s));
            }
            if (SinGamma > 1) SinGamma = 1;
            if (SinGamma < -1) SinGamma = -1;
            if (CosGamma > 1) CosGamma = 1;
            if (CosGamma < -1) CosGamma = -1;

            if (Alfa >= 0)
            {
                Gamma = Math.Asin(SinGamma);
            }
            else
            {
                if (Alfa < 0 && Betta >= 0.0)
                {
                    Gamma = Math.Acos(CosGamma);
                }
                else // Alfa < 0.0 and Betta < 0.0
                {
                    Gamma = -Math.Acos(CosGamma);
                }
            }

            Gamma = Math.Atan2(Betta, Alfa);

            Gamma = _Conv.ToDeg(Gamma);

            if (Gamma<0)
            {
                Gamma = Gamma + 360;
            }
           
            return Gamma;
        }


        // Calculation of rudder angle depending on the angle of rotation around the aircraft axis 
        private Double GetGammaAct(Int32 Number)
        {
            Double GammaAct = 0;
            // By steering wheel number 
            switch (Number)
            {
                case 1:
                    GammaAct = _GammaCalc;
                    break;
                case 2:
                    GammaAct = _GammaCalc + 90;
                    break;
                case 3:
                    GammaAct = _GammaCalc + 180;
                    break;
                case 4:
                    GammaAct = _GammaCalc + 270;
                    break;
            }
            if (_SchemaType)
            {
                GammaAct += 45;
            }
            if (GammaAct > 360)
            {
                GammaAct = GammaAct - 360;
            }
            if (GammaAct < 0)
            {
                GammaAct = GammaAct + 360;
            }
            return GammaAct;
        }

        private Double GetCyT(Double M, Double[] Deltas)
        {

            // When the Mach number exceeds the boundaries of the table, the boundary value is retained 
            if (M < _Cy0Tabl.M.Min())
            {
                M = _Cy0Tabl.M.Min();
            }
            if (M > _Cy0Tabl.M.Max())
            {
                M = _Cy0Tabl.M.Max();
            }

            _Cy0 = _Interp.Interp2D(_Cy0Tabl.M, _Cy0Tabl.Alfa_s, _Cy0Tabl.Data, M, _Alfa_s);

            _dCyd1 = _Interp.Interp4D(_dCyTabl.GammaCalc, _dCyTabl.M, _dCyTabl.Alfa_s, _dCyTabl.Delta, _dCyTabl.Data, _GammaActs[0], M, _Alfa_s, Deltas[0]) ;
            _dCyd2 = _Interp.Interp4D(_dCyTabl.GammaCalc, _dCyTabl.M, _dCyTabl.Alfa_s, _dCyTabl.Delta, _dCyTabl.Data, _GammaActs[1], M, _Alfa_s, Deltas[1]) ;
            _dCyd3 = _Interp.Interp4D(_dCyTabl.GammaCalc, _dCyTabl.M, _dCyTabl.Alfa_s, _dCyTabl.Delta, _dCyTabl.Data, _GammaActs[2], M, _Alfa_s, Deltas[2]);
            _dCyd4 = _Interp.Interp4D(_dCyTabl.GammaCalc, _dCyTabl.M, _dCyTabl.Alfa_s, _dCyTabl.Delta, _dCyTabl.Data, _GammaActs[3], M, _Alfa_s, Deltas[3]) ;

            _CyT = _Cy0 + _dCyd1 + _dCyd2 + _dCyd3 + _dCyd4;

            return _CyT;
        }
        private Double GetCzT(Double M, Double[] Deltas)
        {

            //  When the Mach number exceeds the boundaries of the table, the boundary value is retained 
            if (M < _Cz0Tabl.M.Min())
            {
                M = _Cz0Tabl.M.Min();
            }
            if (M > _Cz0Tabl.M.Max())
            {
                M = _Cz0Tabl.M.Max();
            }
            _Cz0 = _Interp.Interp3D(_Cz0Tabl.GammaCalc, _Cz0Tabl.M, _Cz0Tabl.Alfa_s, _Cz0Tabl.Data, _GammaCalc, M, _Alfa_s);

            _dCzd1 = _Interp.Interp4D(_dCzTabl.GammaCalc, _dCzTabl.M, _dCzTabl.Alfa_s, _dCzTabl.Delta, _dCzTabl.Data, _GammaActs[0], M, _Alfa_s, Deltas[0]);
            _dCzd2 = _Interp.Interp4D(_dCzTabl.GammaCalc, _dCzTabl.M, _dCzTabl.Alfa_s, _dCzTabl.Delta, _dCzTabl.Data, _GammaActs[1], M, _Alfa_s, Deltas[1]);
            _dCzd3 = _Interp.Interp4D(_dCzTabl.GammaCalc, _dCzTabl.M, _dCzTabl.Alfa_s, _dCzTabl.Delta, _dCzTabl.Data, _GammaActs[2], M, _Alfa_s, Deltas[2]);
            _dCzd4 = _Interp.Interp4D(_dCzTabl.GammaCalc, _dCzTabl.M, _dCzTabl.Alfa_s, _dCzTabl.Delta, _dCzTabl.Data, _GammaActs[3], M, _Alfa_s, Deltas[3]);

            _CzT = _Cz0 + _dCzd1 + _dCzd2 + _dCzd3 + _dCzd4;
            return _CzT;
        }
       

        private Double GetmyT(Double M, Double[] Deltas)
        {

            //  When the Mach number exceeds the boundaries of the table, the boundary value is retained 
            if (M < _my0Tabl.M.Min())
            {
                M = _my0Tabl.M.Min();
            }
            if (M > _my0Tabl.M.Max())
            {
                M = _my0Tabl.M.Max();
            }
            _my0 = _Interp.Interp3D(_my0Tabl.GammaCalc, _my0Tabl.M, _my0Tabl.Alfa_s, _my0Tabl.Data, _GammaCalc, M, _Alfa_s);

            _dmyd1 = _Interp.Interp4D(_dmyTabl.GammaCalc, _dmyTabl.M, _dmyTabl.Alfa_s, _dmyTabl.Delta, _dmyTabl.Data, _GammaActs[0], M, _Alfa_s, Deltas[0]);
            _dmyd2 = _Interp.Interp4D(_dmyTabl.GammaCalc, _dmyTabl.M, _dmyTabl.Alfa_s, _dmyTabl.Delta, _dmyTabl.Data, _GammaActs[1], M, _Alfa_s, Deltas[1]);
            _dmyd3 = _Interp.Interp4D(_dmyTabl.GammaCalc, _dmyTabl.M, _dmyTabl.Alfa_s, _dmyTabl.Delta, _dmyTabl.Data, _GammaActs[2], M, _Alfa_s, Deltas[2]);
            _dmyd4 = _Interp.Interp4D(_dmyTabl.GammaCalc, _dmyTabl.M, _dmyTabl.Alfa_s, _dmyTabl.Delta, _dmyTabl.Data, _GammaActs[3], M, _Alfa_s, Deltas[3]);

            _myT = _my0 + _dmyd1 + _dmyd2 + _dmyd3 + _dmyd4;
            return _myT;
        }
        private Double GetmzT(Double M, Double[] Deltas)
        {


            //  When the Mach number exceeds the boundaries of the table, the boundary value is retained 
            if (M < _mz0Tabl.M.Min())
            {
                M = _mz0Tabl.M.Min();
            }
            if (M > _mz0Tabl.M.Max())
            {
                M = _mz0Tabl.M.Max();
            }
            _mz0 = _Interp.Interp2D(_mz0Tabl.M, _mz0Tabl.Alfa_s, _mz0Tabl.Data, M, _Alfa_s);

            _dmzd1 = _Interp.Interp4D(_dmzTabl.GammaCalc, _dmzTabl.M, _dmzTabl.Alfa_s, _dmzTabl.Delta, _dmzTabl.Data, _GammaActs[0], M, _Alfa_s, Deltas[0]);
            _dmzd2 = _Interp.Interp4D(_dmzTabl.GammaCalc, _dmzTabl.M, _dmzTabl.Alfa_s, _dmzTabl.Delta, _dmzTabl.Data, _GammaActs[1], M, _Alfa_s, Deltas[1]);
            _dmzd3 = _Interp.Interp4D(_dmzTabl.GammaCalc, _dmzTabl.M, _dmzTabl.Alfa_s, _dmzTabl.Delta, _dmzTabl.Data, _GammaActs[2], M, _Alfa_s, Deltas[2]);
            _dmzd4 = _Interp.Interp4D(_dmzTabl.GammaCalc, _dmzTabl.M, _dmzTabl.Alfa_s, _dmzTabl.Delta, _dmzTabl.Data, _GammaActs[3], M, _Alfa_s, Deltas[3]);

            _mzT = _mz0 + _dmzd1 + _dmzd2 + _dmzd3 + _dmzd4;

            return _mzT;
        }
        private Double GetmxT(Double M, Double[] Deltas)
        {


            //  When the Mach number exceeds the boundaries of the table, the boundary value is retained 
            if (M < _mx0Tabl.M.Min())
            {
                M = _mx0Tabl.M.Min();
            }
            if (M > _mx0Tabl.M.Max())
            {
                M = _mx0Tabl.M.Max();
            }
            _mx0 = _Interp.Interp3D(_mx0Tabl.GammaCalc, _mx0Tabl.M, _mx0Tabl.Alfa_s, _mx0Tabl.Data, _GammaCalc, M, _Alfa_s);

            
            _dmxd1 = _Interp.Interp4D(_dmxTabl.GammaCalc, _dmxTabl.M, _dmxTabl.Alfa_s, _dmxTabl.Delta, _dmxTabl.Data, _GammaActs[0], M, _Alfa_s, Deltas[0]) ;
            _dmxd2 = _Interp.Interp4D(_dmxTabl.GammaCalc, _dmxTabl.M, _dmxTabl.Alfa_s, _dmxTabl.Delta, _dmxTabl.Data, _GammaActs[1], M, _Alfa_s, Deltas[1]) ;
            _dmxd3 = _Interp.Interp4D(_dmxTabl.GammaCalc, _dmxTabl.M, _dmxTabl.Alfa_s, _dmxTabl.Delta, _dmxTabl.Data, _GammaActs[2], M, _Alfa_s, Deltas[2]) ;
            _dmxd4 = _Interp.Interp4D(_dmxTabl.GammaCalc, _dmxTabl.M, _dmxTabl.Alfa_s, _dmxTabl.Delta, _dmxTabl.Data, _GammaActs[3], M, _Alfa_s, Deltas[3]);

            _mxT = _mx0 + _dmxd1 + _dmxd2 + _dmxd3 + _dmxd4;
            return _mxT;
        }
      
    }
    // Atmosphere
    class Atmosphere
    {
        // Atmospheric tables. Temperature, pressure and air density values as a function of altitude.  
        // Altitude  
        private Double[] _H_tabl = {0, 1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 10e3, 12e3, 14e3, 16e3, 18e3, 20e3,
                      22e3, 24e3, 26e3, 28e3, 30e3, 32e3, 34e3, 36e3, 38e3, 40e3, 42e3, 44e3, 46e3, 48e3, 
	                  50e3, 52e3, 54e3, 56e3, 58e3, 60e3, 62e3, 64e3, 66e3, 68e3, 70e3, 72e3, 74e3, 76e3,
                      78e3, 80e3};
        // Temperature
        private Double[] _T_tabl = {299.660, 293.665, 287.682, 283.656, 276.979, 270.304, 263.632, 256.961, 250.292, 
                      243.626, 236.961, 223.639, 210.325, 197.019, 198.779, 206.713, 214.641, 218.858,
                      222.817, 226.774, 230.728, 236.092, 241.621, 247.147 ,252.670, 258.188, 262.728,
                      267.059, 271.387, 272.350, 272.350, 271.254, 266.544, 261.009, 255.129, 249.253,
                      242.758, 235.906, 229.063, 222.786, 216.928, 211.074, 205.223, 203.225, 201.278,
                      199.331};
        // Pressure 
        private Double[] _P_tabl = {101325.0, 90328.69, 80338.49, 71320.21, 63163.01, 55775.44, 49100.97, 43085.99, 
                      37679.74, 32834.18, 28503.98, 21221.04, 15518.27, 11120.46, 7864.157, 5629.715,
                      4082.017, 2987.766, 2199.570, 1628.370, 1212.014, 907.5227, 684.2101, 519.2440,
                      396.5297, 304.6371, 235.2931, 182.5355, 142.2088, 111.1163, 86.84371, 67.87593, 
                      52.89756, 41.03663, 31.65824, 24.27992, 18.50269, 13.99430, 10.49963, 7.812449,
                      5.767794, 4.223868, 3.066854, 2.213526, 1.592946, 1.142926};
        // Air density 
        private Double[] _ro_tabl = {1.177987, 1.071548, 0.9728581, 0.8759075, 0.7944263, 0.718833, 0.6488298, 0.5841269,
                       0.5244434, 0.469506, 0.4190503, 0.3305654, 0.2570343, 0.1966313, 0.1378225, 0.09487626,
                       0.06625209, 0.04755795, 0.0343896, 0.02501488, 0.01829974, 0.01339103, 0.009864886,
                       0.007319035, 0.00546715, 0.004110401, 0.003119903, 0.002381105, 0.001825475,
                       0.001421309, 0.001110834, 0.000871722, 0.00069136, 0.000547714, 0.000432279,0.000339347,
                       0.000265527, 0.000206657, 0.000159682, 0.000122163, 9.26259E-05, 6.97131E-05, 5.20601E-05,
	                   3.79441E-05, 2.75704E-05, 1.99747E-05, 
	                  };
        // Speed of sound
        private Double[] _a_tabl;

        // Current value of the Mach number 
        private Double _M;
        // Current air density value 
        private Double _ro;
        // Current sound speed value
        private Double _a;

        private Interpolation _Interp;

        public Atmosphere()
        {
            Int16 i;
            _a_tabl = new Double[_H_tabl.Length];
            for (i = 0; i < _H_tabl.Length; i++)
            {
                _a_tabl[i] = 20.0463 * Math.Sqrt(_T_tabl[i]);
            }
            _Interp = new Interpolation(); 

        }
        //Input parameters - altitude and speed
        public Double GetM(Double H, Double V)
        {
            _a = _Interp.Interp1D(_H_tabl, _a_tabl, H);
            _M = V / _a;
            return _M;
        }
        public Double Get_ro(Double H)
        {
            _ro = _Interp.Interp1D(_H_tabl, _ro_tabl, H);
            return _ro;
        }
    }



    // Inertia characteristics
    class InertialFeatures 
    {
        
        private Double[] _m_tabl = {800.0,2161.0,2684.6,3695.8,3838.8,3900.0};
        private Double[] _Jx_tabl = { 3, 130,170,220,230,234};
        private Double[] _Jy_tabl = { 35, 6850, 15060, 22470, 25070, 26060 };
        private Double[] _Jz_tabl = { 35, 6850, 15060, 22470, 25070, 26060 };
             
        private Double _Jx;
        private Double _Jy;
        private Double _Jz;


        private Interpolation _Interp;

        
       public InertialFeatures()
       {
          _Interp = new Interpolation();           
       }

       public Double Get_Jx(Double m)
       {
           _Jx = _Interp.Interp1D(_m_tabl, _Jx_tabl, m);
           return _Jx;
       }

        public Double Get_Jy(Double m)
        {
            _Jy = _Interp.Interp1D(_m_tabl, _Jy_tabl, m);
            return _Jy;
        }
        public Double Get_Jz(Double m)
        {
            _Jz = _Interp.Interp1D(_m_tabl, _Jz_tabl, m);
            return _Jz;
        }


    }

    // Centring characteristics
    class CentringFeatures
    {
        
        private Double[] _m_tabl = { 800.0, 2161.0, 2684.6, 3695.8, 3838.8, 3900.0 };
        private Double[] _Xc_tabl = { 1.727, 3.285, 4.061, 4.447, 4.686, 4.741 };

        private Double _Xc;
        private Double _dXc;
        // Position of centre of mass at the initial moment 
        private const Double Xc0 = 1.727;
      
        private Interpolation _Interp;

        public Double Xc
        {
            get
            {
                return _Xc;
            }
        }

        public CentringFeatures()
        {
            _Interp = new Interpolation();
        }
        private Double Get_Xc(Double m)
        {
            _Xc = _Interp.Interp1D(_m_tabl, _Xc_tabl, m);
            return _Xc;
        }
        public Double Get_dXc(Double m)
        {
            _Xc = Get_Xc(m);
            _dXc = _Xc - Xc0;
            return _dXc;
        }

    }

    // Bivariate table with arrays of values of independent variables 
    class FileData2D
    {
        // Values for the first argument 
        Double[] _M;
        // The values of the second argument 
        Double[] _Alfa_s;
        // Data
        Double[,] _Data;

        public Double[] M
        {
            get
            {
                return _M;
            }
    
        }
        public Double[] Alfa_s
        {
            get
            {
                return _Alfa_s;
            }
            
        }
        public Double[,] Data
        {
            get
            {
                return _Data;
            }

        }

        public FileData2D(Double[,] Data, Double[] M, Double[] Alfa_s)
        {
           
            _Data = Data; 
            _M = M;
            _Alfa_s = Alfa_s; 
 
        }
        

    }
    // Three-dimensional array with arrays of independent variable values 
    class FileData3D
    {
        // Values for the first argument 
        Double[] _GammaCalc;
        //Values for the second argument 
        Double[] _M;
        // Values for the third argument 
        Double[] _Alfa_s;        
        // Data
        Double[,,] _Data;

        public Double[] M
        {
            get
            {
                return _M;
            }

        }
        public Double[] GammaCalc
        {
            get
            {
                return _GammaCalc;
            }

        }
        
        public Double[] Alfa_s
        {
            get
            {
                return _Alfa_s;
            }

        }
        public Double[, ,] Data
        {
            get
            {
                return _Data;
            }

        }
        public FileData3D(Double[,,] Data, Double[] GammaCalc, Double[] M, Double[] Alfa_s)
        {

            _Data = Data;
            _GammaCalc = GammaCalc;
            _M = M;
            _Alfa_s = Alfa_s;

        }


    } 
    // Четырехмерный массив с массивами значений независимых переменных 
    class FileData4D
    {
        //  Values for the first argument 
        Double[] _GammaCalc;
        //  Values for the second argument 
        Double[] _M;
        //  Values for the third argument 
        Double[] _Alfa_s;
        //  Values for the fourth argument 
        Double[] _Delta;

        // Data
        Double[,,,] _Data;

        public Double[] M
        {
            get
            {
                return _M;
            }

        }
        public Double[] GammaCalc
        {
            get
            {
                return _GammaCalc;
            }

        }
        public Double[] Delta
        {
            get
            {
                return _Delta ;
            }

        }
        public Double[] Alfa_s
        {
            get
            {
                return _Alfa_s;
            }

        }
        public Double[,,,] Data
        {
            get
            {
                return _Data;
            }

        }


        public FileData4D(Double[,, ,] Data, Double[] GammaCalc, Double[] M, Double[] Alfa_s, Double[] Delta)
        {

            _Data = Data;
            _GammaCalc = GammaCalc;
            _M = M;
            _Alfa_s = Alfa_s;
            _Delta = Delta;

        }


    }

    // Reading data from the tables 
    class ReadData
    {
        // Checking whether a character is valid for a number 
        private Boolean IsOKforNumber(Char Character)
        {
            Boolean Flag = false;
            if ((Character >= '0' && Character <= '9') || Character == '.' || Character == '-')
            {
                Flag = true;
            }
            return Flag;
        }
        // Converting a data string into an array of numbers 
        private Double[] StringToDouble(String Buf)
        {
            Int16 i,j, WordsCounter;
            Double[] Data_r = new Double[20];
            Double[] Data;  
            String WordS=""; 
            Int16[] WordsLengths=new Int16[20]; 
            j=0;
            i=0;
            
            WordsCounter = 0;
            j = 0;
            Buf = Buf.Trim();
             
            while (true)
                {
                    if (!(IsOKforNumber(Buf[i])) && !(IsOKforNumber(Buf[i+1]))) 
                    {
                        i++;
                        continue;
                    }
                    // Digit
                    if (IsOKforNumber(Buf[i]))
                    {
                        
                        WordS = WordS.Insert(j,Convert.ToString(Buf[i]));   
                        j++;

                    }
                    // Separator or end of line       
                    if (!(IsOKforNumber(Buf[i])) || i == Buf.Length - 1) 
                        
                    {
                            Data_r[WordsCounter] = Convert.ToDouble(WordS);
                            // A new word is recorded 
                            // New word length 
                            WordsLengths[WordsCounter] = j;
                            // New words counter  
                            WordsCounter++;
                            // The index for the next word entry is set to zero 
                            j = 0;
                            WordS = "";
                            if (i == Buf.Length - 1)
                            {
                                break;
                            }
                         
                        
                    }
                    
                   i++;
                    // End of line 
                  
                }
            Data = new Double[WordsCounter];
            for (i = 0; i < WordsCounter; i++)
                Data[i] = Data_r[i];
            return Data;
                   
        }
        // Reading data into a structure with a two-dimensional data table 
        public FileData2D ReadDataFrom2DTable(String FileName)
        {
            StreamReader ReadFile;
            String Buf;
            FileData2D Data2D;
            Double [,] Data; 
            Double[] Var1, Var2, Line; 
            Int16 n,m,i,j;

            ReadFile = new StreamReader(FileName);
            // First headline  
            Buf = ReadFile.ReadLine(); // counted the line
            // Array of values for the first argument 
            Var1 = StringToDouble(ReadFile.ReadLine());
            // Second headline  
            Buf = ReadFile.ReadLine(); // считали линию
            // Array of values for the second argument 
            Var2 = StringToDouble(ReadFile.ReadLine());
            // Number of rows
            n = Convert.ToInt16(Var1.Length);
            // Number of columns
            m = Convert.ToInt16(Var2.Length);
            Data = new Double[n, m];
            // Third headline  
            Buf = ReadFile.ReadLine();
            for (i=0;i<n;i++)
            {
                Line = StringToDouble(ReadFile.ReadLine());
                for (j=0;j<m;j++)
                 {
                   Data[i,j] = Line[j];                                                              
                 }
            }
            Data2D = new FileData2D( Data, Var1, Var2);            
            return Data2D;
  
        }
        // Reading data into a structure with a three-dimensional data table 
        public FileData3D ReadDataFrom3DTable(String FileName)
        {
            StreamReader ReadFile;
            String Buf;
            FileData3D Data3D;
            Double[,,] Data;
            Double[] Var1, Var2, Var3, Line;
            Int16 n, m, k, i, j,p;

            ReadFile = new StreamReader(FileName);
            // First headline  
            Buf = ReadFile.ReadLine(); // считали линию
            // Array of values for the first argument 
            Var1 = StringToDouble(ReadFile.ReadLine());
            // Second heading  
            Buf = ReadFile.ReadLine(); // считали линию
            // Array of values for the second argument 
            Var2 = StringToDouble(ReadFile.ReadLine());
            // Third headline  
            Buf = ReadFile.ReadLine(); // считали линию
            // Array of values for the third argument 
            Var3 = StringToDouble(ReadFile.ReadLine());

            // Height (third dimension)
            k = Convert.ToInt16(Var1.Length);
            // Number of rows
            n = Convert.ToInt16(Var2.Length);
            // Number of columns
            m = Convert.ToInt16(Var3.Length);

            // Three-dimensional data array 
            Data = new Double[k, n, m];

            // Fourth headline  
            Buf = ReadFile.ReadLine();
            for (p = 0; p < k; p++)
            {
                for (i = 0; i < n; i++)
                {
                    Line = StringToDouble(ReadFile.ReadLine());
                    for (j = 0; j < m; j++)
                    {
                        Data[p,i, j] = Line[j];
                    }
                }
            }
            Data3D = new FileData3D(Data, Var1, Var2, Var3);

            return Data3D;
        }
        // Reading data into a structure with a four-dimensional data table 
        public FileData4D ReadDataFrom4DTable(String FileName)
        {
            StreamReader ReadFile;
            String Buf;
            FileData4D Data4D;
            Double[,,,] Data;
            Double[] Var1, Var2, Var3, Var4, Line;
            Int16 n, m, k,u, i, j, p, t;

            ReadFile = new StreamReader(FileName);
            // First headline  
            Buf = ReadFile.ReadLine(); // counted the line
            // Array of values for the first argument 
            Var1 = StringToDouble(ReadFile.ReadLine());
            // Second heading  
            Buf = ReadFile.ReadLine(); //  counted the line
            // Array of values for the second argument 
            Var2 = StringToDouble(ReadFile.ReadLine());
            // Third headline  
            Buf = ReadFile.ReadLine(); // считали линию
            // Array of values for the third argument 
            Var3 = StringToDouble(ReadFile.ReadLine());
            // Fourth headline  
            Buf = ReadFile.ReadLine(); // считали линию
            // Array of values for the fourth argument 
            Var4 = StringToDouble(ReadFile.ReadLine());

            // Fourth dimension
            u =  Convert.ToInt16(Var1.Length);
            // Height (third dimension)
            k = Convert.ToInt16(Var2.Length);
            // Number of rows
            n = Convert.ToInt16(Var3.Length);
            // Number of columns
            m = Convert.ToInt16(Var4.Length);

            // four-dimensional data array 
            Data = new Double[u, k, n, m];

            // Fifth headline  
            Buf = ReadFile.ReadLine();
            for (t = 0; t < u; t++)
            {
                for (p = 0; p < k; p++)
                {
                    for (i = 0; i < n; i++)
                    {
                        Line = StringToDouble(ReadFile.ReadLine());
                        for (j = 0; j < m; j++)
                        {
                            Data[t, p, i, j] = Line[j];
                        }
                    }
                }
            }
            Data4D = new FileData4D(Data, Var1, Var2, Var3,Var4);

            return Data4D;
        }
 
    }

   
}


