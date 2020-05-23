#include "../external/ALGLIB/stdafx.h"
#include "../external/ALGLIB/optimization.h"
// #include "../ V3F.hpp"
#include "../Math/V3F.h"

// time optimal path planner
// using alglib::real_1d_array;
// using alglib::real_2d_array;

// real_1d_array time_optimal_path_planner(V3F my_pos, V3F speed_w, V3F omega, V3F tgt_pos, V3F tgt_angle)
V3F time_optimal_path_planner(V3F my_pos, V3F speed_w, V3F omega, V3F tgt_pos, V3F tgt_angle)
{
  clock_t startTime = clock();
  double g = 9.80665;
  // orientation of the target
  float yaw = tgt_angle.x;
  double ts = sin(yaw); // sin
  double tc = cos(yaw); // cos

  // Cost matrix
  alglib::real_2d_array H =
        "[[2,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,2,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,2,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"
        "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]";

  alglib::real_1d_array f = "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]";
  alglib::real_1d_array s = "[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]";

  alglib::integer_1d_array ct = "[0,0,0,0,0,0,0,0,0,-1,0,0]"; // 0 is =; -1 is <=; 1 is >=
  // constraints
  alglib::real_2d_array c;

  alglib::real_1d_array x;
  alglib::minqpstate state;
  alglib::minqpreport rep;

  // create solver, set quadratic/linear terms
  alglib::minqpcreate(15, state);
  alglib::minqpsetquadraticterm(state, H);
  alglib::minqpsetlinearterm(state, f);
  alglib::minqpsetscale(state, s);
  alglib::minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);

  alglib::real_1d_array x0 = "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]";
  alglib::minqpsetstartingpoint(state, x0);

  for(int i = 0; i < 32; i = i + 1)
  { 
    double t = exp(-2.5 + 0.175 * i);

    // construct the constraints
    double _c[] = {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,my_pos.x, // x of start point
                   0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,my_pos.y, // y of start point
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,my_pos.z, // z of start point
                   pow(t,4),pow(t,3),pow(t,2),pow(t,1),1,0,0,0,0,0,0,0,0,0,0,tgt_pos.x, // x of final point
                   0,0,0,0,0,pow(t,4),pow(t,3),pow(t,2),pow(t,1),1,0,0,0,0,0,tgt_pos.y, // y of final point
                   0,0,0,0,0,0,0,0,0,0,pow(t,4),pow(t,3),pow(t,2),pow(t,1),1,tgt_pos.z, // z of final point
                   0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,speed_w.x, // vx of start point
                   0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,speed_w.y, // vy of start point
                   0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,speed_w.z, // vz of start point
                   tc*4*pow(t,3),tc*3*pow(t,2),tc*2*pow(t,1),tc,0,ts*4*pow(t,3),ts*3*pow(t,2),ts*2*pow(t,1),ts,0,0,0,0,0,0,2, // vx of final point (<= 2) can be tunned
                   -ts*4*pow(t,3),-ts*3*pow(t,2),-ts*2*pow(t,1),-ts,0,tc*4*pow(t,3),tc*3*pow(t,2),tc*2*pow(t,1),tc,0,0,0,0,0,0,0, // vy of final point
                    0,0,0,0,0,0,0,0,0,0,4*pow(t,3),3*pow(t,2),2*pow(t,1),1,0,0}; // vz of final point

    c.setcontent(12,16,_c);

    alglib::minqpsetlc(state, c, ct);

    // Solve problem with BLEIC-based QP solver.
    alglib::minqpoptimize(state);
    alglib::minqpresults(state, x, rep);

    // find the most aggressive trajectory
    double pitch = abs(atan2(sqrt(pow(2 * x(2),2) + pow(2 * x(7),2)), abs(2 * x(12) - g)));
    printf("\n#%i\t\t%+0.6f\t\t%+0.6f\t\t%d",i+1,pitch,sqrt(pow(2 * x(2),2) + pow(2 * x(7),2) + pow(2 * x(12) - g,2)),int(rep.terminationtype));     // display for iteration

    // boundary
    bool throttle_ub_sat = ((2 * x(12) - g < 0) && (sqrt(pow(2 * x(2),2) + pow(2 * x(7),2) + pow(2 * x(12) - g,2)) < 2 * g));
    bool throttle_lb_sat = ((2 * x(12) - g > 0) && (sqrt(pow(2 * x(2),2) + pow(2 * x(7),2) + pow(2 * x(12) - g,2)) < g));

    if(pitch <= M_PI/12 && ( throttle_lb_sat || throttle_ub_sat ))
    {
      printf("\n------------------------------------------------------------------------");
      printf("\nTime for arrival:\t%0.6f",t);
      // clock stop
      clock_t endTime = clock();
      clock_t clockTicksTaken = endTime - startTime;
      double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
      printf("\nOptimizationTime:\t%0.6f",timeInSeconds);
      break;
    }
  }

  // return x;
  return V3F(2 * x[2], 2 * x[7], 2 * x[12] - g);
}


// Construct control input
// V4D minimum_snap_control(mymav mav, mymav tgt, alglib::real_1d_array x, double ts) { 
  // gravity drag force
  // double g = 9.80665;
  // desired acceleration
  // double x_pp = 2 * x[2];
  // double y_pp = 2 * x[7];
  // double z_pp = 2 * x[12] - g;
  // double norm = sqrt(pow(x_pp, 2) + pow(y_pp, 2) + pow(z_pp, 2));
  // store throttle
//   input.throttle = -norm;
//   Eigen::Vector3d zb(x_pp,y_pp,z_pp);
//   zb = - zb / norm;
//   // calculate the desired yaw angle
//   double toward = atan2(tgt.position(1)-mav.position(1),tgt.position(0)-mav.position(0));
//   double increment = toward - mav.angle(0);
//   // normalize the increment yaw angle
//   increment = fmod((increment + M_PI),(2 * M_PI)) - M_PI;

//   if (abs(increment) >=  M_PI/(2/ts))   // assume rotational velocity < 90 degree per second
//   {
//     //increment = abs(increment)/increment * 1/2 * M_PI / (1/ts);
//     increment = abs(increment)/increment * M_PI / (2/ts);
//   }

//   double yaw = mav.angle(0) + increment;
//   Eigen::Vector3d xc(cos(yaw), sin(yaw), 0);
//   Eigen::Vector3d yb;

//   yb = zb.cross(xc);

//   // vector normalization
//   norm = sqrt(pow(yb(0),2) + pow(yb(1),2) + pow(yb(2),2));
//   yb = yb / norm;
//   Eigen::Vector3d xb;
//   xb = yb.cross(zb);

//   // generate the euler angle
//   yaw = atan2(xb(1),xb(0));
//   double pitch = atan2(-xb(2),sqrt(pow(yb(2),2)+pow(zb(2),2)));

//   //normalize pitch
//   pitch = fmod((pitch + M_PI),(2 * M_PI)) - M_PI;
//   double roll = atan2(yb(2),zb(2));

//   // normalize roll
//   pitch = fmod((pitch + M_PI),(2 * M_PI)) - M_PI;

//   // store yaw, pitch, roll
//   input.yaw = yaw;
//   input.pitch = pitch;
//   input.roll = roll;

//   // display the result
//   printf("\n\tthrottle:\t%+0.6f",input.throttle);
//   printf("\n\tyaw:\t\t%+0.6f",input.yaw);
//   printf("\n\tpitch:\t\t%+0.6f",input.pitch);
//   printf("\n\troll:\t\t%+0.6f",input.roll);
  // return V4D();
// }
