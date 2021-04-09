using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;


public class rotation : MonoBehaviour
{

    public Vector3 point, aixs;

    private const double rotate_once = 0.5;
    private int status;
    private double fixedDeltaTime;
    private double angle, target, delta, virtual_angle;
    private int times;
    // Start is called before the first frame update
    void Start()
    {
        status = 0;
        angle = target = delta = 0;
        times = 0;
        fixedDeltaTime = (double)UnityEngine.Time.fixedDeltaTime;
     }

    private void set_target(double[] param)
    {
        status = 0;
        target = param[0];
        times = (int) (param[1] / fixedDeltaTime);

        delta = (target - angle) / times; //the tiny-angle every time rotate

        virtual_angle = angle;
        status = (int) param[2];

    }
    // Update is called once per frame
    void Update()
    {
    }
    public int get_status()
    {
        return status;
    }
    void Impuls(int dir)
    {
        gameObject.transform.RotateAround(point, aixs, (float)rotate_once * dir);
        angle += dir * rotate_once;
    }
    private void FixedUpdate()
    {
        if (status == 0) return;
        if (times == 0)
        {
            if (status == 2)
            {
                double[] param = new double[3];
                param[0] = 0;
                param[1] = 0.08;
                param[2] = 1.1;
                set_target(param);
                return;
            }   else
            if (status == 1)
            {
                status = 0;
                return;
            }
            return;
        }
        times--;
        double new_angle = virtual_angle + delta;
        if ((int)(new_angle / rotate_once) != (int)(virtual_angle / rotate_once))
        {
            Impuls(Math.Sign(delta));
        }
        virtual_angle = new_angle;


    }
}
