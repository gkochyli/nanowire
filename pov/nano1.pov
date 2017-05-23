#include "colors.inc"
 
#declare substr12 = 
pigment {
	color rgb<0,1,0>}
 
#declare substr3 = 
pigment {
	color rgb<0,0,1>}
	 
#declare nanowire = 
pigment {
	color rgb<1,1,1>}


camera {
   location <110,8,-8>
   look_at <50,0,50>
   // angle 10 // direction <0, 0,  1.7>
   //right x*image_width/image_height
   //look_at <20,0.5,2>
   //location<-5,3,0>
   //look_at <0,0,0>
}

background{ color rgb <0,0,0> }

//light_source {< 150, 10, -50> White }
light_source {< 150, 10, 30> White }
//#declare vecCupos = <0.0,0.0,0.0>;
#fopen sub1 "1_2layer.txt" read
#while (defined(sub1))
    #declare vecCupos = <0.0,0.0,0.0>;
    #read (sub1,vecCupos)
    sphere {vecCupos, 0.5 
        texture {
            pigment {
             substr12
            }
        }
    }
#end
#fclose sub1


#fopen sub3 "3rd_layer.txt" read
#while (defined(sub3))
    #declare vecCupos = <0.0,0.0,0.0>;
    #read (sub3,vecCupos)
    sphere {vecCupos, 0.5 
        texture {
            pigment {
             substr3
            }
        }
    }
#end
#fclose sub3


#fopen nano "nanowire_coordinates.txt" read
#while (defined(nano))
	#declare vecCupos = <0.0,0.0,0.0>;
	#read (nano,vecCupos)
	sphere {vecCupos, 0.66
		texture {
			pigment { nanowire }
			}
		}
#end
#fclose nano