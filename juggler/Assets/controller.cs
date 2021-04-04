using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics;
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
    public GameObject[] engine;
    public double dev;
    public double P, I, D;
    public double hP;

    private const int Buffsize = 100;
    private double Hit_Height = 8;
    private double Plate_Length = 299;
    private double Middle_Length = 80;
    private double Lower_Length = 89;
    private double engine_pos = 80.477;

    private double fixedDeltaTime;
    private myRandom ND = new myRandom();
    private int index;
    private double[] h, x, y;

    private double predict_time_delta = 0.03;
    private double predict_time_min = 0.01;
    private double predict_time_max = 0.03;
    

    // Start is called before the first frame update
    void Start()
    {
        h = new double[Buffsize];
        x = new double[Buffsize];
        y = new double[Buffsize];
        for (int i = 0; i < Buffsize; i++)
            h[i] = -1;
        index = -1;
        fixedDeltaTime = Time.fixedDeltaTime;
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    double Clip(double a, double l, double r)
    {
        a = Math.Min(a, r);
        a = Math.Max(a, l);
        return a;
    }

    private double calc_phi(double Height, double phi)
    {
        Vector2 p1 = new Vector2((float) (-Math.Cos(phi) * Plate_Length / 2), (float)(-Math.Sin(phi) * Plate_Length / 2 + Height));
        Vector2 p2 = new Vector2((float) -engine_pos, 0);
        float R1 = (float)Middle_Length;
        float R2 = (float)Lower_Length;
        
        return 0;
    }
    private void FixedUpdate()
    {
        int status = 0;
        for (int i = 0; i < 4; i++)
            status = Math.Max(status, engine[i].GetComponent<rotation>().get_status());
        if (status == 2) return;


        index = (index + 1) % Buffsize;
        h[index] = ball.transform.position.y;
        x[index] = ball.transform.position.x;
        y[index] = ball.transform.position.z - ball.transform.localScale.x / 2;

        
        
        h[index] += ND.Next_nd(0, dev);
        x[index] += ND.Next_nd(0, dev);
        y[index] += ND.Next_nd(0, dev);

        if (status == 1) return;

        int cnt = 0;
        for (int i = 0; i < Buffsize; i++)
        {
            int new_index = (index - i + Buffsize) % Buffsize;
            if (h[new_index] == -1) break;
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

        double[] res_x = Fit.Polynomial(new_t, new_x, 1);
        double[] res_y = Fit.Polynomial(new_t, new_y, 1);
        double[] res_h = Fit.Polynomial(new_t, new_h, 2);

        double count = cnt - 1 + predict_time_delta / fixedDeltaTime;

        double h_predict = count * count * res_h[2] + count * res_h[1] + res_h[0];
        if (h_predict < Hit_Height)
        { 
            double c1 = (-res_h[1] + Math.Sqrt(res_h[1] * res_h[1] - 4 * res_h[2] * (res_h[0] - Hit_Height))) / 2 / res_h[2];
            double c2 = (-res_h[1] - Math.Sqrt(res_h[1] * res_h[1] - 4 * res_h[2] * (res_h[0] - Hit_Height))) / 2 / res_h[2];
            if (Math.Abs(c1 - count) < Math.Abs(c2 - count))
                count = c1;
            else
                count = c2;

            if (double.IsInfinity(count) || double.IsNaN(count)) count = cnt - 1 + predict_time_delta / fixedDeltaTime;
                else count = Clip(count, cnt - 1 + predict_time_min / fixedDeltaTime, cnt - 1 + predict_time_max / fixedDeltaTime);

            double pre_x = count * res_x[1] + res_x[0];
            double pre_y = count * res_y[1] + res_y[0];

            double pre_vh = 2 * res_h[2] * count + res_h[1];
            double pre_vx = res_x[1];
            double pre_vy = res_y[1];

            double phi_x = (Math.PI / 2 - Math.Atan2(-pre_vh, pre_vx)) / 2;
            double phi_y = (Math.PI / 2 - Math.Atan2(-pre_vh, pre_vy)) / 2;


            double[] param = new double[3];
            param[0] = calc(Hit_Height, phi_x);
            param[1] = 0.03;
            param[2] = 2.1;
            engine[2].SendMessage("set_target", param);

        }


        for (int i = 0; i < 4; i++) {
            double[] param = new double[3];
            param[0] = 7;
            param[1] = 0.05;
            param[2] = 2.1;
            engine[i].SendMessage("set_target", param);
        }

    }
}
