      subroutine get_cmd_args(ierr)
c
c     $Id$
c
c**********************************************************************
c     Purpose: Interrogate the command line for specific options.
c
c     Usage:
c
c     call get_cmd_args (ierr)
c
c     ierr     : (out) 0 if no errors, 1 if some error occurred.
c
c     Original coding: Thomas R. Faulkner, MRJ Technology Solutions
c**********************************************************************

      character*8 arg
c
      common /proces/ numprocs
c
c     Begin get_cmd_args.
c
c     Initialize common data.
c
      numprocs = -1
c
c     Get the total number of command line arguments.
c
      num_args = iargc()
      if (num_args .eq. 2) then
#if defined CRAY_TIME
         call pxfgetarg (1, arg, 3, ierror)
#else
         call getarg (1, arg)
#endif
         if (arg .eq. "-np" .or. arg .eq. "-n") then
#if defined CRAY_TIME
            call pxfgetarg (2, arg, 4, ierror)
#else
            call getarg (2, arg)
#endif
            read (arg, 100, err = 110) n
 100        format (i3)
            numprocs = n
            ierr = 0
            go to 111
 110        continue
            ierr = 1
 111        continue
         else
            ierr = 1
         end if
      else
         ierr = 1
      end if
      return
      end
