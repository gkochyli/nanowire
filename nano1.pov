#include "colors.inc"
 
#declare substr12 = 
pigment {
     color rgb<0,1,0>}
 
#declare substr3 = 
pigment {
     color rgb<0,0,1>}
 
camera {
   //location <10, -12, 0>   
   location <50,-20,0.5>
   //angle 45 // direction <0, 0,  1.7>
   //right x*image_width/image_height
   //look_at <20,0.5,2> 
   look_at <50,50,100>
}                      

light_source {< 0, 12, -30> White }
//light_source {< 10, 10, 10> White }
#declare vecCupos = <0.0,0.0,0.0>;
#fopen sub1 "1st_layer.txt" read
#while (defined(sub1))
    #read (sub1,vecCupos)
    sphere {vecCupos, 1 
        texture {
            pigment {
             substr12
            }
        }
    }
#end
