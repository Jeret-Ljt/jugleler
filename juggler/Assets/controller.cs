using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Double;

using System.IO;
using System.Text;

using System;

public class myRandom
{
    private System.Random rd = new System.Random();
    private bool V_hot = false;
    private double V, U;
    public double Next_nd(double mean, double dev)
    {
        double Up;
        if (V_hot)
        {
            V_hot = false;
            Up = V;
        }
        else
        {
            double S;
            do
            {
                U = (rd.NextDouble() - 0.5) / 0.5;
                V = (rd.NextDouble() - 0.5) / 0.5;
                S = U * U + V * V;
            } while (S >= 1 || S == 0);
            double Fp = Math.Sqrt(-2 * Math.Log(S) / S);
            V = V * Fp;
            V_hot = true;
            U = U * Fp;
            Up = U;
        }
        return Up * dev + mean;
    }
}


public class controller : MonoBehaviour
{
    public GameObject ball;
    public GameObject Cube;
    public double dev = 0;
    public double P = 0.3, I = 0, D = 0;
    public double Pz = 0.1;
    public double Target_H = 50, Targe_H_vz = 0;
    

    private const double CubeDepth = 0.01;
    private const double BallRadius = 1.5;
    private const int Buffsize = 10000;
    private const double predict_time_delta = 0.05;

    private myRandom ND = new myRandom();
    private double fixedDeltaTime;
    private int index;
    private double[] h, x, y;
    private double top_h;
    private double move_up_more_time;
    private int status;
    private double vz, x_phi, y_phi;
    private int collide_tick = -100;

    private double vx = 0, vy = 0;

    private Vector3 eulerAngle_in_collide;
    private System.Random rd = new System.Random();


    private StreamWriter filexy;
    private StreamWriter fileh;

    int collide_times = 0;

    
    void Clear()
    {
        top_h = 0;
        index = -1;
        status = -1;
        move_up_more_time = 0.05;
        for (int i = 0; i < Buffsize; i++) h[i] = -100;   
    }

    
    void OnCollisionEnter(Collision collision)
    {
        if (collide_times == 10000) return;
        collide_times++;
        if (status == 1)
        {
            filexy.WriteLine(ball.transform.position.x);
            filexy.WriteLine(ball.transform.position.z);
            fileh.WriteLine(top_h - Target_H);
            fileh.Flush();
            filexy.Flush();
            // Debug.Log(123);
            Debug.Log(top_h);
            Debug.Log(ball.transform.position.x);
            Debug.Log(ball.transform.position.z);
        }
        else
        if (status == -1)
        {
            if (collide_tick >= 0)
            {
                int cnt = 0;
                for (int i = 0; i < Buffsize; i++)
                {
                    int new_index = (index - i + Buffsize) % Buffsize;
                    if (h[new_index] == -100) break;
                    cnt++;
                }

//                Debug.Log(cnt);
//                Debug.Log(collide_tick);
                double[] new_h = new double[collide_tick + 1];
                double[] new_x = new double[collide_tick + 1];
                double[] new_y = new double[collide_tick + 1];
                double[] new_t = new double[collide_tick + 1];


                for (int i = 0; i <= collide_tick; i++)
                {
                    int original_index = i;
                    int new_index = i;
                    new_h[new_index] = h[original_index];
                    new_x[new_index] = x[original_index];

                   // Debug.Log(new_x[new_index]);
                    new_y[new_index] = y[original_index];
                    new_t[new_index] = new_index * fixedDeltaTime;
                }

                double[] res_x = Fit.Polynomial(new_t, new_x, 1);
                double[] res_y = Fit.Polynomial(new_t, new_y, 1);
                double[] res_h = Fit.Polynomial(new_t, new_h, 2);

                Vector3 v0 = new Vector3((float)res_x[1], (float)(res_h[2] * 2 * (collide_tick) * fixedDeltaTime + res_h[1]), (float)res_y[1]);

                double[] new_h2 = new double[cnt - collide_tick];
                double[] new_x2 = new double[cnt - collide_tick];
                double[] new_y2 = new double[cnt - collide_tick];
                double[] new_t2 = new double[cnt - collide_tick];
                for (int i = collide_tick; i < cnt; i++)
                {
                    int original_index = i;
                    int new_index = i - collide_tick;
                    new_h2[new_index] = h[original_index];
                    new_x2[new_index] = x[original_index];


                //    Debug.Log(new_x2[new_index]);
                    new_y2[new_index] = y[original_index];
                    new_t2[new_index] = new_index * fixedDeltaTime;
                }


                double[] res_x2 = Fit.Polynomial(new_t2, new_x2, 1);
                double[] res_y2 = Fit.Polynomial(new_t2, new_y2, 1);
                double[] res_h2 = Fit.Polynomial(new_t2, new_h2, 2);

                Vector3 v1 = new Vector3((float)res_x2[1], (float)(res_h2[1]), (float)res_y2[1]);


               
                double cosy = Math.Cos(eulerAngle_in_collide.z / 180 * Math.PI);
                double siny = Math.Sin(eulerAngle_in_collide.z / 180 * Math.PI);

                double cosx = Math.Cos(eulerAngle_in_collide.x / 180 * Math.PI);
                double sinx = Math.Sin(eulerAngle_in_collide.x / 180 * Math.PI);

                double[,] R = {
                            { cosy, sinx*siny, -cosx*siny },
                            { 0, cosx,sinx },
                            {siny,-sinx*cosy, cosx*cosy }
                          };
                var mR = DenseMatrix.OfArray(R);
              
                double[,] V0 =
                {
                {v0.x - Cube.GetComponent<Rigidbody>().velocity.x},
                {v0.z - Cube.GetComponent<Rigidbody>().velocity.z},
                {v0.y - Cube.GetComponent<Rigidbody>().velocity.y},
            };
                double[,] V1 =
                {
                {v1.x - Cube.GetComponent<Rigidbody>().velocity.x},
                {v1.z - Cube.GetComponent<Rigidbody>().velocity.z},
                {v1.y - Cube.GetComponent<Rigidbody>().velocity.y},
            };


                var mV0 = DenseMatrix.OfArray(V0);
                var mV1 = DenseMatrix.OfArray(V1);

             
                var mR_V0 = mR * mV0;
                var mR_V1 = mR * mV1;

                Vector3 new_v0 = new Vector3((float)mR_V0.At(0, 0), (float)(mR_V0.At(2, 0) ), (float)mR_V0.At(1, 0));
                Vector3 new_v1 = new Vector3((float)mR_V1.At(0, 0), (float)(mR_V1.At(2, 0) ), (float)mR_V1.At(1, 0));

                Vector3 v2 = new Vector3(new_v1.x - new_v0.x * 0.98f , new_v1.y - new_v0.y * (-0.81f), new_v1.z - new_v0.z * 0.98f);

                Vector3 ratiao = new Vector3(new_v1.x / new_v0.x, new_v1.y / new_v0.y, new_v1.z / new_v0.z);

                Debug.Log(ratiao);

                //  Debug.Log(v2);
                //  fileh.WriteLine(v2.y.ToString("0.000"));
                 fileh.WriteLine(new_v0.y.ToString("0.00") + " " + new_v1.y.ToString("0.00"));
                filexy.WriteLine(new_v0.x.ToString("0.00") + " " + new_v1.x.ToString("0.00"));
              //  filexy.WriteLine(v2.x);

                fileh.Flush();
                filexy.Flush();

               
                //   Debug.Log(new_v0);
               //    Debug.Log(new_v1);
                 Debug.Log(ratiao.y.ToString("0.000"));



                collide_tick = -100;
                ball.transform.position = new Vector3(0, (float)(rd.NextDouble() * Target_H), 0);
                ball.GetComponent<Rigidbody>().velocity = new Vector3((float)(rd.NextDouble() - 0.5) * 300, 0, (float)(rd.NextDouble() - 0.5) * 300);

              
        
                Targe_H_vz = rd.NextDouble() * 100;
                vx = rd.NextDouble() * 100;
                vy = rd.NextDouble() * 100;

                Cube.transform.position = new Vector3(0, -10, 0);
                Cube.transform.localEulerAngles = new Vector3((float)rd.NextDouble() * 45, 0, (float)rd.NextDouble() * 45);
                Cube.GetComponent<Rigidbody>().velocity = new Vector3((float)vx, (float)Targe_H_vz, (float)vy);
                Cube.GetComponent<Rigidbody>().angularVelocity = new Vector3(0, 0, 0);

              //  ball.GetComponent<Rigidbody>().angularVelocity = new Vector3((float)(rd.NextDouble() - 0.5) * 100, (float)(rd.NextDouble() - 0.5) * 100, (float)(rd.NextDouble() - 0.5) * 100);
                

                return;
            }
            else
            {
                if (collide_tick == -100)
                {
                    
                    collide_tick = -99;
                }
                else
                if (collide_tick == -99)
                {
                    index = -1;
                    for (int i = 0; i < Buffsize; i++) h[i] = -100;
                    collide_tick = -98;
                }   else
                if (collide_tick == -98)
                {
                    collide_tick = index;
                    //                Debug.Log(h[index]); 
                    eulerAngle_in_collide = -Cube.transform.rotation.eulerAngles;
                }
            }
        }
    }


    // Start is called before the first frame update
    void Start()
    {
        filexy = new StreamWriter("./xydata.txt");
        fileh = new StreamWriter("./hdata.txt");

        h = new double[Buffsize];
        x = new double[Buffsize];
        y = new double[Buffsize];
        fixedDeltaTime = Time.fixedDeltaTime;
        move_up_more_time = 0.1;

        for (int i = 0; i < Buffsize; i++) h[i] = -100;
        index = -1;
        top_h = 0;
        status = -1; //0 : predict status, 1: move status,  -1: test status


        
        ball.transform.position = new Vector3(0, (float)Target_H, 0);
        ball.GetComponent<Rigidbody>().velocity = new Vector3((float)(rd.NextDouble() - 0.5) * 300, 0, (float)(rd.NextDouble() - 0.5) * 300);
        ball.GetComponent<Rigidbody>().maxAngularVelocity = 3;
      //  ball.GetComponent<Rigidbody>().angularVelocity = new Vector3(10, 10, 10);
      //  Cube.transform.position = new Vector3(0, -10, 0);
      //  Cube.transform.rotation = new Quaternion(0, 0, 0, 0);
       // Cube.GetComponent<Rigidbody>().velocity = new Vector3(0, 0, 0);
       // Cube.GetComponent<Rigidbody>().angularVelocity = new Vector3(0, 0, 0);

    }

    // Update is called once per frame
    void Update()
    {
        
    }



    private void FixedUpdate()
    {
        if (collide_times == 10000)
        {
            filexy.Close();
            fileh.Close();
            return;
        }
            if (status == 0)
        {
            index = (index + 1) % Buffsize;
            h[index] = ball.transform.position.y;
            x[index] = ball.transform.position.x;
            y[index] = ball.transform.position.z;
            top_h = Math.Max(top_h, h[index]);

            h[index] += ND.Next_nd(0, dev);
            x[index] += ND.Next_nd(0, dev);
            y[index] += ND.Next_nd(0, dev);
            
            int cnt = 0;
            for (int i = 0; i < Buffsize; i++)
            {
                int new_index = (index - i + Buffsize) % Buffsize;
                if (h[new_index] == -100) break;
                cnt++;
            }

            double[] new_h = new double[cnt];
            double[] new_x = new double[cnt];
            double[] new_y = new double[cnt];
            double[] new_t = new double[cnt];
            for (int i = 0; i < cnt; i++)
            {
                int original_index = (index - i + Buffsize) % Buffsize;
                int new_index = cnt - 1 - i;
                new_h[new_index] = h[original_index];
                new_x[new_index] = x[original_index];
                new_y[new_index] = y[original_index];
                new_t[new_index] = new_index;
            }

            if (cnt <= 5) return;
            double[] res_x = Fit.Polynomial(new_t, new_x, 1);
            double[] res_y = Fit.Polynomial(new_t, new_y, 1);
            double[] res_h = Fit.Polynomial(new_t, new_h, 2);

            double count = cnt - 1 + predict_time_delta / fixedDeltaTime; //预计碰撞时间
            double h_predict = count * count * res_h[2] + count * res_h[1] + res_h[0];


            if (h_predict - BallRadius <= 0)
            {
                double t = (-res_h[1] - Math.Sqrt(res_h[1] * res_h[1] - 4 * res_h[2] * (res_h[0] - BallRadius))) / (2 * res_h[2]);

                count = t;
                double new_predict_time_delta = (count - (cnt - 1)) * fixedDeltaTime;

                Debug.Log(new_predict_time_delta);
                double x_predict = count * res_x[1] + res_x[0];
                double y_predict = count * res_y[1] + res_y[0];

                y_phi = x_predict * P + res_x[1] / fixedDeltaTime * D;
                x_phi = -(y_predict * P + res_y[1] / fixedDeltaTime * D);

                vz = (Target_H - top_h) * Pz + Targe_H_vz;
                Debug.Log(vz);

                Cube.transform.localEulerAngles = new Vector3((float)x_phi, 0, (float)y_phi);
                //unit = cm/s

                Cube.transform.position = new Vector3((float)x_predict, (float)(-vz * new_predict_time_delta - CubeDepth / 2), (float)y_predict);
                Cube.GetComponent<Rigidbody>().velocity = new Vector3(0, (float)vz, 0);
                status = 1; // moving plate
            }
        }
        else if (status == 1)
        {
            int last_index = index;
            index = (index + 1) % Buffsize;
            h[index] = ball.transform.position.y;
            x[index] = ball.transform.position.x;
            y[index] = ball.transform.position.z;


            if (h[index] > h[last_index])
            {
               // Debug.Log(top_h);
                //Debug.Log(h[last_index]);
                move_up_more_time -= fixedDeltaTime;
                if (move_up_more_time <= 0)
                {
                    Clear();
                    Cube.transform.position = new Vector3(0, -10, 0);
                    Cube.GetComponent<Rigidbody>().velocity = new Vector3(0, 0, 0);
                    Cube.GetComponent<Rigidbody>().angularVelocity = new Vector3(0, 0, 0);
                    Cube.transform.localEulerAngles = new Vector3(0, 0, 0);
                }
            }
        }
        else if (status == -1)
        {
            
            index = (index + 1) % Buffsize;
            h[index] = ball.transform.position.y;
            x[index] = ball.transform.position.x;
            y[index] = ball.transform.position.z;

            
     
        }
    }
}
