# Make sure that you are building openblas with OPENBLAS_DYNAMIC_ARCH=1
# You may have to wipe your openblas build to ensure that it is built
# with support for all architectures, or else performance may suffer.

JULIAHOME := $(abspath ../../..)
include ../../../Make.inc

JULIA_PKGDIR:=$(shell echo $$JULIA_PKGDIR)
ifeq ($(JULIA_PKGDIR),)
	JULIA_PKGDIR:=$(shell echo ~)/.julia
endif

DMG_VERSION_SUFFIX:=$(shell [ $$(git describe --tags --exact-match 2>/dev/null) ] && echo $(JULIA_VERSION) || echo $(JULIA_VERSION)-$(JULIA_COMMIT))

APP_VERSION_SUFFIX:=$(shell echo $(JULIA_VERSION) | grep -o '^[0-9]\+.[0-9]\+')

DMG_NAME:=Julia-$(DMG_VERSION_SUFFIX).dmg
APP_NAME:=Julia-$(APP_VERSION_SUFFIX).app
VOL_NAME:=Julia

all: clean dmg

dmg:
	make -C ../../../deps install-git
	make -C ../../.. binary-dist
	tar zxf ../../../julia-*.tar.gz
	mv julia-* julia
	-mkdir -p ./julia/libexec ./julia/share
	-cp -a $(build_libexecdir)/git* ./julia/libexec
	-cp -a $(build_datarootdir)/git* ./julia/share
	rm -f julia/lib/*.{a,la}
	-mkdir dmg
	platypus -a Julia -p /bin/bash -V $(JULIA_VERSION) -R -u "The Julia Project" -i julia.icns -Q julia.icns -o "None" -I org.julialang -x -f julia script ./dmg/$(APP_NAME)
	-codesign -s "AFB379C0B4CBD9DB9A762797FC2AB5460A2B0DBE" --deep ./dmg/Julia-$(APP_VERSION_SUFFIX).app
	-cp -f julia.icns dmg/.VolumeIcon.icns
	-mkdir dmg/.background
	-cp -f background.pdf dmg/.background/
	-ln -fs /Applications ./dmg/Applications
	-chmod 755 ./dmg/$(APP_NAME)/Contents/MacOS/Julia
	hdiutil create temp.dmg -size 500m -ov -volname "$(VOL_NAME)" -srcfolder dmg -format UDRW
	-mkdir mnt
	hdiutil attach -noverify -readwrite -noautoopen temp.dmg -mountpoint mnt
#	SetFile -a C mnt
	-osascript window.applescript $(APP_NAME)
	@echo "We're going to chown the .app file to root:admin now, which requires sudo.  You may be asked for your password:"
	-sudo chown -R root:admin mnt/$(APP_NAME)
	sync
	sync
	hdiutil detach mnt
	hdiutil convert temp.dmg -format UDZO -imagekey zlib-level=9 -o "$(DMG_NAME)"


clean:
	@echo "We have to use sudo here to clean out folders owned by root.  You may be asked for your password"
	sudo rm -fr julia dmg *.dmg mnt
