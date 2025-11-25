#include "rigid_body_integrator.hpp"

#include "rigid_body.hpp"

namespace physsim
{
    void explicitEuler(RigidBody& body, double stepSize)
    {
        if(body.type() == RigidBody::EType::Static)
        {
            return;
        }
        // get current position and rotation of the body
        Eigen::Vector3d x    = body.position();
        Eigen::Quaterniond q = body.rotation();

        // get current linear and angular velocity of the body
        Eigen::Vector3d v = body.linearVelocity();
        Eigen::Vector3d w = body.angularVelocity();

        // update position of the body using the linear velocity and update body accordingly
        Eigen::Vector3d xnew = x + stepSize * v;
        body.setPosition(xnew);

        // quaternion-based angular velocity update of rotation and update body accordingly
        Eigen::Quaterniond wq(0, w.x(), w.y(), w.z());
        Eigen::Quaterniond qnew = (q + 0.5 * stepSize * wq * q).normalized();
        body.setRotation(qnew);

        // get current linear and angular momentum of body
        Eigen::Vector3d p = body.linearMomentum();
        Eigen::Vector3d l = body.angularMomentum();

        // get force and torque that are currently applied to body
        Eigen::Vector3d f = body.force();
        Eigen::Vector3d t = body.torque();

        // compute new linear momentum
        Eigen::Vector3d pnew = p + stepSize * f;

        // convert from linear momentum to linear velocity and update the body accordingly
        Eigen::Vector3d vnew = body.massInverse() * pnew;
        body.setLinearVelocity(vnew);

        // compute new angular momentum
        Eigen::Matrix3d I    = body.inertiaWorld();
        Eigen::Vector3d lnew = l + stepSize * (t - w.cross(I * w));

        // convert from angular momentum to angular velocity and update the body accordingly
        Eigen::Vector3d wnew = body.inertiaWorldInverse() * lnew;
        body.setAngularVelocity(wnew);
    }

    void symplecticEuler(RigidBody& body, double stepSize)
    {
        
        if(body.type() == RigidBody::EType::Static)
        {
            return;
        }
        // TODO: get current position and rotation of the body
        Eigen::Vector3d x = body.position();
        Eigen::Quaterniond q = body.rotation();

        // TODO: get current linear and angular momentum of body
        Eigen::Vector3d p = body.linearMomentum();
        Eigen::Vector3d l = body.angularMomentum();

        // TODO: get force and torque that are currently applied to body
        Eigen::Vector3d f = body.force();
        Eigen::Vector3d t = body.torque();

        // TODO: compute new linear momentum
        Eigen::Vector3d p_new = p + stepSize * f;

    
        // TODO: convert from linear momentum to linear velocity and update the body accordingly
        Eigen::Vector3d v_new = p_new * body.massInverse();
        body.setLinearVelocity(v_new);

        // TODO: compute new angular momentum
        Eigen::Vector3d w     = body.angularVelocity();
        Eigen::Matrix3d iner  = body.inertiaWorld();
        Eigen::Vector3d l_new = l + stepSize * (t - w.cross(iner * w));


        // TODO: convert from angular momentum to angular velocity and update the body accordingly
        Eigen::Vector3d w_new = body.inertiaWorldInverse() * l_new;
        body.setAngularVelocity(w_new);

        // TODO: update position of the body using the linear velocity and update body accordingly
        Eigen::Vector3d x_new = x + stepSize * v_new;
        body.setPosition(x_new);

        // TODO: quaternion-based angular velocity update of rotation and update body accordingly
        Eigen::Quaterniond wq(0, w_new.x(), w_new.y(), w_new.z());
        Eigen::Quaterniond q_new = (q + 0.5 * stepSize * wq * q).normalized();
        body.setRotation(q_new);

        // Eigen::FullPivHouseholderQR() // use for qr decomp
    }


    void implicitEuler(RigidBody& body, double stepSize)
    {
        
        if(body.type() == RigidBody::EType::Static)
        {
            return;
        }
        // See for additional explanations: https://www.gdcvault.com/play/1022196/Physics-for-Game-Programmers-Numerical

        // TODO: get current position and rotation of the body
        Eigen::Vector3d x = body.position();
        Eigen::Quaterniond q = body.rotation();

        // TODO: get current linear and angular momentum of body
        Eigen::Vector3d p = body.linearMomentum();
        Eigen::Vector3d l = body.angularMomentum();

        // TODO: get force and torque that are currently applied to body
        Eigen::Vector3d f = body.force();
        Eigen::Vector3d t = body.torque();

        // TODO: compute new linear momentum
        Eigen::Vector3d p_new = p + stepSize * f;

        // TODO: convert from linear momentum to linear velocity and update the body accordingly
        Eigen::Vector3d v_new = body.massInverse() * p_new;
        body.setLinearVelocity(v_new);

        // TODO: Convert current angular velocity to body coordinates (initial guess wb0)
        Eigen::Vector3d w   = body.angularVelocity();
        Eigen::Vector3d w_b0 = q.inverse() * w;

        // TODO: Compute residual vector f(wb0) from the from angular velocity in body-coordinates
        Eigen::Matrix3d bi = body.inertiaBody();
        Eigen::Vector3d f_wb0 = stepSize * w_b0.cross(bi * w_b0);

        // TODO: Compute the Jacobian of f at wb.
        ///skew in place of omega gen
        Eigen::Matrix3d omega_ib_wb0 = skew(bi * w_b0);
        Eigen::Matrix3d omega_wb0 = skew(w_b0);
        Eigen::Matrix3d j = bi + stepSize * (omega_wb0 * bi - omega_ib_wb0);

        // TODO: Linearly solve for the update step delta_wb, for example using a QR decomposition
        //use fun 
        Eigen::Vector3d delta_wb = j.colPivHouseholderQr().solve(-f_wb0);

        // TODO: Apply the Newton-Raphson iteration by adding delta_wb to the current angular velocity
        Eigen::Vector3d w_b = w_b0 + delta_wb ;        

        // TODO: Transform the angular velocity back to world coordinates
        Eigen::Vector3d w_new = q * w_b;

        // TODO: explicitly integrate the torque and update the body accordingly
        w_new += body.inertiaWorldInverse() * stepSize * t;
        body.setAngularVelocity(w_new);

        // TODO: update position of the body using the linear velocity and update body accordingly
        Eigen::Vector3d x_new = x + stepSize * v_new;
        body.setPosition(x_new);

        // TODO: quaternion-based angular velocity update of rotation
        Eigen::Quaterniond wq(0, w_new.x(), w_new.y(), w_new.z());
        Eigen::Quaterniond q_new = (q + 0.5 * stepSize * wq * q).normalized();
        body.setRotation(q_new);
    }

}
