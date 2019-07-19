{-# LANGUAGE CPP #-}
{-# OPTIONS_GHC -fno-warn-missing-import-lists #-}
{-# OPTIONS_GHC -fno-warn-implicit-prelude #-}
module Paths_primerclip (
    version,
    getBinDir, getLibDir, getDynLibDir, getDataDir, getLibexecDir,
    getDataFileName, getSysconfDir
  ) where

import qualified Control.Exception as Exception
import Data.Version (Version(..))
import System.Environment (getEnv)
import Prelude

#if defined(VERSION_base)

#if MIN_VERSION_base(4,0,0)
catchIO :: IO a -> (Exception.IOException -> IO a) -> IO a
#else
catchIO :: IO a -> (Exception.Exception -> IO a) -> IO a
#endif

#else
catchIO :: IO a -> (Exception.IOException -> IO a) -> IO a
#endif
catchIO = Exception.catch

version :: Version
version = Version [0,3,9,0] []
bindir, libdir, dynlibdir, datadir, libexecdir, sysconfdir :: FilePath

bindir     = "/home/irish/github/primerclip/.stack-work/install/x86_64-linux/lts-11.0/8.2.2/bin"
libdir     = "/home/irish/github/primerclip/.stack-work/install/x86_64-linux/lts-11.0/8.2.2/lib/x86_64-linux-ghc-8.2.2/primerclip-0.3.9.0-GWHGlFCMArZ9UpP4hK9nN6-primerclip-test"
dynlibdir  = "/home/irish/github/primerclip/.stack-work/install/x86_64-linux/lts-11.0/8.2.2/lib/x86_64-linux-ghc-8.2.2"
datadir    = "/home/irish/github/primerclip/.stack-work/install/x86_64-linux/lts-11.0/8.2.2/share/x86_64-linux-ghc-8.2.2/primerclip-0.3.9.0"
libexecdir = "/home/irish/github/primerclip/.stack-work/install/x86_64-linux/lts-11.0/8.2.2/libexec/x86_64-linux-ghc-8.2.2/primerclip-0.3.9.0"
sysconfdir = "/home/irish/github/primerclip/.stack-work/install/x86_64-linux/lts-11.0/8.2.2/etc"

getBinDir, getLibDir, getDynLibDir, getDataDir, getLibexecDir, getSysconfDir :: IO FilePath
getBinDir = catchIO (getEnv "primerclip_bindir") (\_ -> return bindir)
getLibDir = catchIO (getEnv "primerclip_libdir") (\_ -> return libdir)
getDynLibDir = catchIO (getEnv "primerclip_dynlibdir") (\_ -> return dynlibdir)
getDataDir = catchIO (getEnv "primerclip_datadir") (\_ -> return datadir)
getLibexecDir = catchIO (getEnv "primerclip_libexecdir") (\_ -> return libexecdir)
getSysconfDir = catchIO (getEnv "primerclip_sysconfdir") (\_ -> return sysconfdir)

getDataFileName :: FilePath -> IO FilePath
getDataFileName name = do
  dir <- getDataDir
  return (dir ++ "/" ++ name)
