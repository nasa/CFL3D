!  ---------------------------------------------------------------------------
!  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
!  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
!  with point-matched, patched, overset, or embedded connectivities. Both
!  multigrid and mesh sequencing are available in time-accurate or
!  steady-state modes.
!
!  Copyright 2001 United States Government as represented by the Administrator
!  of the National Aeronautics and Space Administration. All Rights Reserved.
! 
!  The CFL3D platform is licensed under the Apache License, Version 2.0 
!  (the "License"); you may not use this file except in compliance with the 
!  License. You may obtain a copy of the License at 
!  http://www.apache.org/licenses/LICENSE-2.0. 
! 
!  Unless required by applicable law or agreed to in writing, software 
!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!  License for the specific language governing permissions and limitations 
!  under the License.
!  ---------------------------------------------------------------------------
!
!
!   Reynolds Stress Model Periodic BC Module
!
!   called by: bc2005.F, bc2005i_d.F, bc2005j_d.F, bc2005k_d.F
!
!   Main functions:
!      1) construct the rotation matrix for tensor rotation
!      2) apply the boundary conditions for 7 turbulence variables
!            a) 6-stress components
!            b) 1-length scale equation component (omega, epsilon, or zeta ...)

MODULE module_stm_2005
  IMPLICIT NONE
  private :: tensor_rotate
CONTAINS
  ! build the rotation martix and its transpose
  ! based on the rotation angle specified, c.f. the 
  ! notes in the bc2005.F
  SUBROUTINE stm2k5_get_rotmat(thx, thy, thz, rn,rnt)
    implicit none
    real, intent(in) :: thx, thy, thz
    real, intent(out) :: rn(3,3),rnt(3,3)

    if(abs(thx)>0) then
       rn(:,1) = (/1., 0., 0./)
       rn(:,2) = (/0., cos(thx), sin(thx)/)
       rn(:,3) = (/0.,-sin(thx), cos(thx)/)
    else if(abs(thy)>0) then
       rn(:,1) = (/cos(thy), 0. , -sin(thy)/)
       rn(:,2) = (/0.,  1.0,  0./)
       rn(:,3) = (/sin(thy), 0. ,  cos(thy)/)
    else
       rn(:,1) = (/ cos(thz), sin(thz), 0./)
       rn(:,2) = (/-sin(thz), cos(thz), 0./)
       rn(:,3) = (/0., 0., 1.0/)
    endif
    ! get the transpose of the rotation matrix
    rnt(:,1) = rn(1,:)
    rnt(:,2) = rn(2,:)
    rnt(:,3) = rn(3,:)

    return
  END SUBROUTINE stm2k5_get_rotmat

  ! apply the boundary condition
  SUBROUTINE stm2k5_bc(vin, rn, rnt, vout)
    REAL, INTENT(in) :: vin(7),rn(3,3),rnt(3,3)
    REAL, intent(out) :: vout(7)

    CALL tensor_rotate(vin,rn,rnt,vout)
    vout(7) = vin(7)
  END SUBROUTINE stm2k5_bc

  ! rotate the tensor.
  ! Tout = R *  Tin  * R_t  
  ! which is not the same as that in the symmetry bc
  ! Here we rotate the tensor itself, not the reference frame.
  
  SUBROUTINE tensor_rotate(tin, rn, rnt, tout)
    implicit none
    REAL, INTENT(in) :: tin(6),rn(3,3),rnt(3,3)
    real, intent(out) :: tout(6)

    real :: t33(3,3),s(3,3),v(3,3)

    t33(1,1)= tin(1)
    t33(2,2)= tin(2)
    t33(3,3)= tin(3)

    t33(1,2)= tin(4)
    t33(2,3)= tin(5)
    t33(1,3)= tin(6)

    t33(2,1)= t33(1,2)
    t33(3,2)= t33(2,3)
    t33(3,1)= t33(1,3)

    s(:,1) = t33(:,1)*rnt(1,1) + t33(:,2)*rnt(2,1) + t33(:,3)*rnt(3,1)
    s(:,2) = t33(:,1)*rnt(1,2) + t33(:,2)*rnt(2,2) + t33(:,3)*rnt(3,2)
    s(:,3) = t33(:,1)*rnt(1,3) + t33(:,2)*rnt(2,3) + t33(:,3)*rnt(3,3)

    v(:,1) = rn(:,1)*s(1,1) + rn(:,2)*s(2,1) + rn(:,3)*s(3,1)
    v(:,2) = rn(:,1)*s(1,2) + rn(:,2)*s(2,2) + rn(:,3)*s(3,2)
    v(:,3) = rn(:,1)*s(1,3) + rn(:,2)*s(2,3) + rn(:,3)*s(3,3)

    tout(1) = v(1,1)
    tout(2) = v(2,2)
    tout(3) = v(3,3)
    tout(4) = v(1,2)
    tout(5) = v(2,3)
    tout(6) = v(1,3)
    return
  end subroutine tensor_rotate
END MODULE module_stm_2005
