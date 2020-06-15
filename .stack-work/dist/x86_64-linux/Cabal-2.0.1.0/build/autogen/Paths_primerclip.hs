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
version = Version [0,3,9,1] []
bindir, libdir, dynlibdir, datadir, libexecdir, sysconfdir :: FilePath

bindir     = "/home/admin/stack/primerclip/.stack-work/install/x86_64-linux/3f0c863dd7d182fb30b0b70be660352fc230b53ac8c4b8bf581e074e77625bd8/8.2.2/bin"
libdir     = "/home/admin/stack/primerclip/.stack-work/install/x86_64-linux/3f0c863dd7d182fb30b0b70be660352fc230b53ac8c4b8bf581e074e77625bd8/8.2.2/lib/x86_64-linux-ghc-8.2.2/primerclip-0.3.9.1-ENS9d9mgpiN1AkwRgsXwFu"
dynlibdir  = "/home/admin/stack/primerclip/.stack-work/install/x86_64-linux/3f0c863dd7d182fb30b0b70be660352fc230b53ac8c4b8bf581e074e77625bd8/8.2.2/lib/x86_64-linux-ghc-8.2.2"
datadir    = "/home/admin/stack/primerclip/.stack-work/install/x86_64-linux/3f0c863dd7d182fb30b0b70be660352fc230b53ac8c4b8bf581e074e77625bd8/8.2.2/share/x86_64-linux-ghc-8.2.2/primerclip-0.3.9.1"
libexecdir = "/home/admin/stack/primerclip/.stack-work/install/x86_64-linux/3f0c863dd7d182fb30b0b70be660352fc230b53ac8c4b8bf581e074e77625bd8/8.2.2/libexec/x86_64-linux-ghc-8.2.2/primerclip-0.3.9.1"
sysconfdir = "/home/admin/stack/primerclip/.stack-work/install/x86_64-linux/3f0c863dd7d182fb30b0b70be660352fc230b53ac8c4b8bf581e074e77625bd8/8.2.2/etc"

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
