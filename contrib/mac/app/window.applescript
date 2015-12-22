#!/usr/bin/osascript

on run argv
   set appname to item 1 of argv

   tell application "Finder"
      set currpath to container of (path to me) as alias
      set mnt to folder "mnt" of currpath

      tell mnt
    	 open
	 tell container window
	    set current view to icon view
            set toolbar visible to false
	    set statusbar visible to false
            set bounds to {400, 100, 700, 400}
	 end tell

	 set viewopts to the icon view options of container window
	 set arrangement of viewopts to not arranged
	 set icon size of viewopts to 72	    
	 set background picture of viewopts to file ".background:background.pdf"

	 set position of item appname of container window to {80, 170}
	 set position of item "Applications" of container window to {248, 170}

         update without registering applications

         close
      end tell
   end tell
end run

