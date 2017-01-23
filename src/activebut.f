        SUBROUTINE ACTIVEBUT(JUST,MODECUT)
        IMPLICIT NONE
        INTEGER JUST,MODECUT
C
        INTEGER NCOLOR1,NCOLOR2,NCOLOR3,NCOLOR4,NCOLOR5
C------------------------------------------------------------------------------
        NCOLOR1=-3-1
        NCOLOR2=-2-1
        NCOLOR3=-4-1
        NCOLOR4=-5-1
        NCOLOR5=-7-1
C------------------------------------------------------------------------------
        CALL BUTTON(2,'[s]ave',NCOLOR1)
        CALL BUTTON(12,'[i]math',NCOLOR1)
C------------------------------------------------------------------------------
        CALL BUTTON(3,'[b]oundary',0)
C------------------------------------------------------------------------------
        CALL BUTTON(5,'[z]oom',NCOLOR2)
        CALL BUTTON(6,'[w]hole',NCOLOR2)
        IF(JUST.EQ.1)THEN
          CALL BUTTON(7,'x[=]y',NCOLOR2)
        ELSE
          CALL BUTTON(7,'  x![=]y',NCOLOR2)
        END IF
        IF(MODECUT.EQ.1)THEN
          CALL BUTTON(15,'single',NCOLOR2)
        ELSEIF(MODECUT.EQ.2)THEN
          CALL BUTTON(15,'region',NCOLOR2)
        ELSE
          CALL BUTTON(15,'all',NCOLOR2)
        END IF
        CALL BUTTON(16,'[x] cut',NCOLOR2)
        CALL BUTTON(17,'[y] cut',NCOLOR2)
C------------------------------------------------------------------------------
        CALL BUTTON(8,'[m]apping',NCOLOR5)
C------------------------------------------------------------------------------
        CALL BUTTON(10,'stack 0',0)
C------------------------------------------------------------------------------
        CALL BUTTON(20,'s[Tt]atistic',0)
C------------------------------------------------------------------------------
        CALL BUTTON(22,'postscript',0)
        CALL BUTTON(23,'RGB *.ppm',0)
C------------------------------------------------------------------------------
        CALL BUTTON(24,'[Cc]entroid',0)
C------------------------------------------------------------------------------
        CALL BUTTON(29,'[.]special',0)
C------------------------------------------------------------------------------
        CALL BUTTON(31,'resize ima',0)
C------------------------------------------------------------------------------
        CALL BUTTON(34,'box[9]oper',0)
        CALL BUTTON(44,'[o]ffsets',0)
C------------------------------------------------------------------------------
        CALL BUTTON(161,'zoom',0)
        CALL BUTTON(161,'zoom',NCOLOR3)
        CALL BUTTON(162,'min[,]max',0)
        CALL BUTTON(162,'min[,]max',NCOLOR3)
        CALL BUTTON(163,'z1[/]z2',0)
        CALL BUTTON(163,'z1[/]z2',NCOLOR3)
        CALL BUTTON(164,'BG[:]FG',0)
        CALL BUTTON(164,'BG[:]FG',NCOLOR3)
C------------------------------------------------------------------------------
        CALL BUTTON(41,'zoom',0)
        CALL BUTTON(41,'zoom',NCOLOR4)
        CALL BUTTON(42,'whole',0)
        CALL BUTTON(42,'whole',NCOLOR4)
        CALL BUTTON(43,'min,max',0)
        CALL BUTTON(43,'min,max',NCOLOR4)
        END
