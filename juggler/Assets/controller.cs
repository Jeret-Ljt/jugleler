using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics;



public class controller : MonoBehaviour
{
    public GameObject ball;
    public GameObject engine0, engine1, engine2, engine3;

    private const int Buffsize = 100;
    private int index;
    private double[] h, x, y;

    // Start is called before the first frame update
    void Start()
    {
        h = new double[Buffsize];
        x = new double[Buffsize];
        y = new double[Buffsize];
        for (int i = 0; i < Buffsize; i++)
            h[i] = -1;
        index = -1;
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    private void FixedUpdate()
    {
        index = (index + 1) % Buffsize;
        h[index] = ball.transform.position.y;
        x[index] = ball.transform.position.x;
        y[index] = ball.transform.position.z;



    }
}
