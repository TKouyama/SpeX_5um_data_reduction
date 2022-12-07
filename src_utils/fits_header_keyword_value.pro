;;
;; Fits headerから指定したキーワードのValueをStr形式で得る関数
;;

function fits_header_keyword_value,header,keyword
  keyword_length = strlen(keyword)

  if keyword_length le 8 then begin
    keyword_pos = where(strtrim(strmid(header,0,8),2) eq keyword)
    if keyword_pos[0] ne -1 then begin
      str_keyword_line = header[keyword_pos]
      equal_pos = strpos(str_keyword_line,'=')
      slash_pos = strpos(str_keyword_line,' /')+1

      if slash_pos gt 8 then begin
        value = strmid(str_keyword_line,equal_pos+1,slash_pos-equal_pos-1)
        value = strtrim(value,2)
      endif else begin
        value = strmid(str_keyword_line,equal_pos+1,30)
        value = strtrim(value,2)
      endelse
    endif else begin
      value = ''
    endelse

  endif else begin

    ;result = strcmp(header, keyword)
    ;keyword_pos = where(strtrim(strmid(header,0,keyword_length),2) eq strtrim(keyword,2))
    keyword_pos = where(strmid(header,0,keyword_length) eq keyword)

    if keyword_pos[0] ne -1 then begin
      str_keyword_line = header[keyword_pos]
      equal_pos = strpos(str_keyword_line,'=')
      slash_pos = strpos(str_keyword_line,' /')+1

      if slash_pos gt 8 then begin
        value = strmid(str_keyword_line,equal_pos+1,slash_pos-equal_pos-1)
        value = strtrim(value,2)
      endif else begin
        value = strmid(str_keyword_line,equal_pos+1,30)
        value = strtrim(value,2)
      endelse
    endif else begin
      value = ''
    endelse
    
  endelse



  return,value[0]
end
