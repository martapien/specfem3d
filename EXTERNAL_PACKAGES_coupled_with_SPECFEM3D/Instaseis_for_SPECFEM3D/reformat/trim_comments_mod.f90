module trim_comments

contains

subroutine trim_string(string)
   character(500), intent(inout) :: string
   character(500) :: string_trim
   integer :: ind_hash
   ind_hash = index(string, '#') - 1
   string_trim = string(3:ind_hash)
   string = trim(string_trim)
end subroutine trim_string

end module trim_comments
