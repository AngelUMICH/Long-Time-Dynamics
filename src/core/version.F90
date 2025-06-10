subroutine version
#ifdef GIT_VERSION_F
  write(6,'(2A)') "Git Hash: ", GIT_VERSION_F
#endif
end subroutine version