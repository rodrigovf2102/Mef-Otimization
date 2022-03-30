# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=Portico3d - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to Portico3d - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "Portico3d - Win32 Release" && "$(CFG)" !=\
 "Portico3d - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "Portico3d.mak" CFG="Portico3d - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Portico3d - Win32 Release" (based on\
 "Win32 (x86) Console Application")
!MESSAGE "Portico3d - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "Portico3d - Win32 Debug"
F90=fl32.exe
RSC=rc.exe

!IF  "$(CFG)" == "Portico3d - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
OUTDIR=.
INTDIR=.

ALL : "$(OUTDIR)\Portico3d.exe"

CLEAN : 
	-@erase ".\Portico3d.exe"
	-@erase ".\Jacobi.obj"
	-@erase ".\DefinicaoParametros.obj"
	-@erase ".\Portico3d.obj"
	-@erase ".\Assemble.obj"
	-@erase ".\Solvers.obj"

# ADD BASE F90 /Ox /c /nologo
# ADD F90 /Ox /c /nologo
F90_PROJ=/Ox /c /nologo 
# ADD BASE RSC /l 0x416 /d "NDEBUG"
# ADD RSC /l 0x416 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/Portico3d.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/Portico3d.pdb" /machine:I386 /out:"$(OUTDIR)/Portico3d.exe" 
LINK32_OBJS= \
	"$(INTDIR)/Jacobi.obj" \
	"$(INTDIR)/DefinicaoParametros.obj" \
	"$(INTDIR)/Portico3d.obj" \
	"$(INTDIR)/Assemble.obj" \
	"$(INTDIR)/Solvers.obj"

"$(OUTDIR)\Portico3d.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "Portico3d - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
OUTDIR=.
INTDIR=.

ALL : "$(OUTDIR)\Portico3d.exe"

CLEAN : 
	-@erase ".\Portico3d.exe"
	-@erase ".\Jacobi.obj"
	-@erase ".\DefinicaoParametros.obj"
	-@erase ".\Portico3d.obj"
	-@erase ".\Assemble.obj"
	-@erase ".\Solvers.obj"
	-@erase ".\Portico3d.ilk"
	-@erase ".\Portico3d.pdb"

# ADD BASE F90 /Zi /c /nologo
# ADD F90 /Zi /c /nologo
F90_PROJ=/Zi /c /nologo /Fd"Portico3d.pdb" 
# ADD BASE RSC /l 0x416 /d "_DEBUG"
# ADD RSC /l 0x416 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/Portico3d.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:yes\
 /pdb:"$(OUTDIR)/Portico3d.pdb" /debug /machine:I386\
 /out:"$(OUTDIR)/Portico3d.exe" 
LINK32_OBJS= \
	"$(INTDIR)/Jacobi.obj" \
	"$(INTDIR)/DefinicaoParametros.obj" \
	"$(INTDIR)/Portico3d.obj" \
	"$(INTDIR)/Assemble.obj" \
	"$(INTDIR)/Solvers.obj"

"$(OUTDIR)\Portico3d.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.for.obj:
   $(F90) $(F90_PROJ) $<  

.f.obj:
   $(F90) $(F90_PROJ) $<  

.f90.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "Portico3d - Win32 Release"
# Name "Portico3d - Win32 Debug"

!IF  "$(CFG)" == "Portico3d - Win32 Release"

!ELSEIF  "$(CFG)" == "Portico3d - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\Portico3d.for

"$(INTDIR)\Portico3d.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Assemble.for

"$(INTDIR)\Assemble.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Jacobi.for

"$(INTDIR)\Jacobi.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Entrada.txt

!IF  "$(CFG)" == "Portico3d - Win32 Release"

!ELSEIF  "$(CFG)" == "Portico3d - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Saida.dat

!IF  "$(CFG)" == "Portico3d - Win32 Release"

!ELSEIF  "$(CFG)" == "Portico3d - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Entradainformacoes.txt

!IF  "$(CFG)" == "Portico3d - Win32 Release"

!ELSEIF  "$(CFG)" == "Portico3d - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\DefinicaoParametros.for

"$(INTDIR)\DefinicaoParametros.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Saidaotimizacao.txt

!IF  "$(CFG)" == "Portico3d - Win32 Release"

!ELSEIF  "$(CFG)" == "Portico3d - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Solvers.for

"$(INTDIR)\Solvers.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
