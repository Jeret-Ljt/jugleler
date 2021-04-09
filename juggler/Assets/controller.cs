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
    public double h_P, h_const;
    public double delta_t;
    public RectTransform UI_pos;
    private const double CubeDepth = 0.2;
    private const int Buffsize = 100;
    private const double Hit_Height = 8.5;
    private const double Plate_Length = 29.9;
    private const double Middle_Length = 8.0;
    private const double Lower_Length = 8.9;
    private const double engine_pos = 8.0477;
    private const double target_height = 35;
    private double fixedDeltaTime;
    private myRandom ND = new myRandom();
    private int index;
    private double[] h, x, y;
    private double top_h;
    

    private const double predict_time_delta = 0.08;
    private const double predict_time_min = 0.03;
    private const double predict_time_max = 0.08;

    private const double angle_limit_max = 30;
    private const double angle_limit_min = -30;

    // Start is called before the first frame update
    void Start()
    {
        h = new double[Buffsize];
        x = new double[Buffsize];
        y = new double[Buffsize];
        fixedDeltaTime = Time.fixedDeltaTime;
        for (int i = 0; i < Buffsize; i++)
            h[i] = -1;
        index = -1;
        top_h = -1;
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
    double Normalize(double a, double o, double l, double r)
    {
        if (double.IsInfinity(a) || double.IsNaN(a)) return o;
            else return Clip(a, l, r);
    }

    private double calc_phi(double Height, double phi)
    {
        if (double.IsInfinity(phi) || double.IsNaN(phi)) return 0;
        if (double.IsInfinity(Height) || double.IsNaN(phi)) return 0;

        Vector2 p1 = new Vector2((float) (-Math.Cos(phi) * Plate_Length / 2), (float)(Math.Sin(phi) * Plate_Length / 2 + Height));
        Vector2 p2 = new Vector2((float) -engine_pos, 0);


        float R1 = (float)Middle_Length;
        float R2 = (float)Lower_Length;

        float d = Vector2.Distance(p1, p2);
        float x = (d * d + R1 * R1 - R2 * R2) / 2 / d;
        float y = Mathf.Sqrt(R1 * R1 - x * x);

        Vector2 O = (p2 - p1) / d * x + p1;
        Vector2 f = Vector2.Perpendicular((p2 - p1) / d);
        Vector2 C1 = O + f * y;
        Vector2 C2 = O - f * y;

        Vector2 Original = p2 - new Vector2((float)Lower_Length, 0);
        if (Vector2.Distance(C1, Original) < Vector2.Distance(C2, Original))
        {
            Original = Original - p2;
            C1 = C1 - p2;

            double sign = C1.x * Original.y - C1.y * Original.x;
            if (double.IsNaN(sign) || double.IsInfinity(sign)) return 0;

            return Math.Sign(sign) * Vector2.Angle(C1, Original);
        }
        else
        {
            Original = Original - p2;
            C2 = C2 - p2;

            double sign = C2.x * Original.y - C2.y * Original.x;
            if (double.IsNaN(sign) || double.IsInfinity(sign)) return 0;
            return Math.Sign(sign) * Vector2.Angle(C2, Original);
        }
    }

    void Clear()
    {
        top_h = -1;
        index = -1;
        for (int i = 0; i < Buffsize; i++) h[i] = -1;
    }
    private void FixedUpdate()
    {
        int status = 0;
        for (int i = 0; i < 4; i++)
            status = Math.Max(status, engine[i].GetComponent<rotation>().get_status());
        if (status == 2) return;


        index = (index + 1) % Buffsize;
        h[index] = ball.transform.position.y - ball.transform.localScale.x / 2;
        x[index] = ball.transform.position.x;
        y[index] = ball.transform.position.z;

        
        h[index] += ND.Next_nd(0, dev);
        x[index] += ND.Next_nd(0, dev);
        y[index] += ND.Next_nd(0, dev);


        top_h = Math.Max(top_h, h[index]);

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

        if (cnt <= 5) return;
        double[] res_x = Fit.Polynomial(new_t, new_x, 1);
        double[] res_y = Fit.Polynomial(new_t, new_y, 1);
        double[] res_h = Fit.Polynomial(new_t, new_h, 2);


        double count = cnt - 1 + predict_time_delta / fixedDeltaTime;
        double h_predict = count * count * res_h[2] + count * res_h[1] + res_h[0];


            //(-b+- sqrt(b^2-4ac)) / (2a)
            double c1 = (-res_h[1] + Math.Sqrt(res_h[1] * res_h[1] - 4 * res_h[2] * (res_h[0] - Hit_Height))) / 2 / res_h[2];
            double c2 = (-res_h[1] - Math.Sqrt(res_h[1] * res_h[1] - 4 * res_h[2] * (res_h[0] - Hit_Height))) / 2 / res_h[2];
            if (Math.Abs(c1 - count) < Math.Abs(c2 - count) && c1 >= cnt - 1)
                count = c1;
            else
                count = c2;

            count = Normalize(count,
                cnt - 1 + predict_time_delta / fixedDeltaTime,
                cnt - 1 + predict_time_min / fixedDeltaTime,
                cnt - 1 + predict_time_max / fixedDeltaTime);

            double pre_x = count * res_x[1] + res_x[0];
            double pre_y = count * res_y[1] + res_y[0];
            UI_pos.localPosition = new Vector3((float)pre_x, -(float)pre_y, 0);


        if (h_predict < Hit_Height && status == 0)
        {
            double pre_vh = 2 * res_h[2] * count + res_h[1];
            double pre_vx = res_x[1];
            double pre_vy = res_y[1];

            double phi_x = (Math.Atan2(-pre_vh, pre_vx) - Math.PI / 2) / 2;
            double phi_y = (Math.Atan2(-pre_vh, pre_vy) - Math.PI / 2) / 2;

            double final_phi_x = I * phi_x + P * (-pre_x);
            double final_phi_y = I * phi_y + P * (-pre_y);

            double delta_h = -CubeDepth / 2;

            delta_h += h_P * (target_height - top_h) + h_const;
            delta_h = Clip(delta_h, -CubeDepth / 2, 5);

            Debug.Log(top_h);

            double[] param = new double[3];
            param[1] = (count - (cnt - 1)) * fixedDeltaTime + delta_t;
            param[2] = 2.1;

            final_phi_x = Normalize(final_phi_x, 0, -Math.PI / 3, Math.PI / 3);
            final_phi_y = Normalize(final_phi_y, 0, -Math.PI / 3, Math.PI / 3);

            param[0] = calc_phi(Hit_Height + delta_h, final_phi_x);
            param[0] = Normalize(param[0], 0, angle_limit_min, angle_limit_max);
 
            engine[2].SendMessage("set_target", param);

            param[0] = calc_phi(Hit_Height + delta_h, -final_phi_x);
            param[0] = Normalize(param[0], 0, angle_limit_min, angle_limit_max);

            engine[0].SendMessage("set_target", param);

            param[0] = calc_phi(Hit_Height + delta_h, final_phi_y);
            param[0] = Normalize(param[0], 0, angle_limit_min, angle_limit_max);

            engine[1].SendMessage("set_target", param);

            param[0] = calc_phi(Hit_Height + delta_h, -final_phi_y);
            param[0] = Normalize(param[0], 0, angle_limit_min, angle_limit_max);
           
            engine[3].SendMessage("set_target", param);

            Clear();
        }


    }
}
