using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ball : MonoBehaviour
{

    public GameObject controller;
    private void Init()
    {
        var rb = gameObject.GetComponent<Rigidbody>();
        rb.velocity = new Vector3(0, 0, 0);
        rb.position = new Vector3(-15, 30, -15);
        rb.AddForce(400 * (Vector3.forward + Vector3.right));
    }

    // Start is called before the first frame update
    void Start()
    {
        Init();
        
    }

    // Update is called once per frame
    void Update()
    {
        
    }

   
    private void FixedUpdate()
    {
        Vector3 pos = gameObject.transform.position;
        

    }
}
