module Lib
    ( module Lib
    ) where

import System.Process
import Control.Applicative
import Options.Applicative hiding (flag)
import Data.Semigroup ((<>))
import Control.Monad
import Control.Monad.IO.Class (liftIO)
import Control.Exception (evaluate)
import qualified Data.Either as U
import Data.Ord (comparing)
import qualified Data.Char as C
import Data.Foldable (toList)
import qualified Data.Map.Strict as M
import Data.List
import qualified Data.ByteString.Char8 as B
import qualified Data.Attoparsec.ByteString.Char8 as A
import qualified Conduit as P
import qualified Data.Conduit.Attoparsec as CA
import qualified Data.Vector as V
import Data.Bits
import Data.Digits
import Data.Ord (comparing, Down)
import Data.Monoid hiding ((<>))
import qualified Data.IntMap.Strict as I
import Data.Maybe
import qualified Data.Set as S
import GHC.Generics (Generic)

{--
    Jonathan Irish
    Post-alignment primer trimming tool v0.3.1

--}

------------------------------------------------------------------------------
--------------------------------  Data Types  --------------------------------
------------------------------------------------------------------------------

-- alignment record containing alignment data for each read
data AlignedRead = AlignedRead { qname :: B.ByteString
                               , flag :: Int
                               , rname :: UChr
                               , pos :: Integer
                               , endpos :: Integer
                               , mapqual :: Integer
                               , cigar :: B.ByteString
                               , trimdcigar :: B.ByteString
                               , cigmap :: CigarMap
                               , trimdcigmap :: CigarMap
                               , rnext :: B.ByteString
                               , pnext :: Integer
                               , tlen :: Integer
                               , refseq :: B.ByteString
                               , basequal :: B.ByteString
                               , optfields :: B.ByteString
                               , strand :: B.ByteString
                               , paired :: Bool
                               , mapped :: Bool
                               , fint :: [BedRecord]
                               , rint :: [BedRecord]
                               , pintflag :: Bool
                               , trimdflag :: Bool
                               , trimdpos :: Integer
                               , trimdendpos :: Integer
                               , headerstrings :: Header
                               , isheader :: Bool
                               , mid :: B.ByteString
                               , xmid :: B.ByteString
                               , tbed :: BedRecord
                               , trimdToZeroLength :: Bool
                               } deriving (Eq, Show, Generic)

instance Ord AlignedRead where
    compare =  comparing rname
            <> comparing pos
            <> comparing mapqual

defaultAlignment = AlignedRead { qname = "NONE"
                               , flag = 0
                               , rname = NONE
                               , pos = 0
                               , endpos = 0
                               , mapqual = 0
                               , cigar = "*"
                               , trimdcigar = "*"
                               , cigmap = [(0,"*")]
                               , trimdcigmap = [(0,"*")]
                               , rnext = "*"
                               , pnext = 0
                               , tlen = 0
                               , refseq = ""
                               , basequal = ""
                               , optfields = ""
                               , strand = ""
                               , paired = False
                               , mapped = False
                               , fint = []
                               , rint = []
                               , pintflag = False
                               , trimdflag = False
                               , trimdpos = 0
                               , trimdendpos = 0
                               , headerstrings = []
                               , isheader = False
                               , mid = ""
                               , xmid = ""
                               , tbed = defaultBed
                               , trimdToZeroLength = False
                               }

-- 180321 organize aligned reads by name, with read1 and read2 organized by
-- primary and secondary Alignments
-- NOTE: this requires the input SAM file to be name sorted, and enables the
-- update of the mate-related flag bits and mate alignment start positions
-- for each trimmed alignment
-- (this should remove most if not all ValidateSamFile errors)
data PairedAln = PairedAln { r1prim :: AlignedRead
                           , r2prim :: AlignedRead
                           , r1secs :: [AlignedRead]
                           , r2secs :: [AlignedRead]
                           } deriving (Eq, Show, Generic)

instance Ord PairedAln where
    compare = comparing r1prim
           <> comparing r2prim

defaultPairedAln = PairedAln defaultAlignment defaultAlignment [] []

-- 180206 add BEDPE primer coordinates input file option
data BEDPE = BEDPE { chr1 :: UChr
                   , start1 :: Integer
                   , end1 :: Integer
                   , chr2 :: UChr
                   , start2 :: Integer
                   , end2 :: Integer
                   , bedpename :: B.ByteString
                   } deriving (Show, Eq, Generic)

instance Ord BEDPE where
    compare =  comparing chr1
            <> comparing start1

defaultBEDPE = BEDPE NONE 0 0 NONE 0 0 "DEFAULTBEDPE"

data BedRecord = BedRecord { bedchr :: UChr
                           , bedstart :: Integer
                           , bedend :: Integer
                           , bedname :: B.ByteString
                           } deriving (Show, Eq, Generic)

instance Ord BedRecord where
    compare =  comparing bedchr
            <> comparing bedstart

defaultBed = BedRecord { bedchr = NONE
                       , bedstart = 0
                       , bedend = 0
                       , bedname = "DEFAULT"
                       }

-- 170926 use master file as source of primer and target coords
data MasterRecord = MasterRecord { mchrom :: UChr
                                 , mtargetstart :: Integer
                                 , mtargetend :: Integer
                                 , mampname :: B.ByteString
                                 , mforstart :: Integer
                                 , mforend :: Integer
                                 , mrevstart :: Integer
                                 , mrevend :: Integer
                                 , mforname :: B.ByteString
                                 , mrevname :: B.ByteString
                                 , mforseq :: B.ByteString
                                 , mrevseq :: B.ByteString
                                 } deriving (Show, Eq)

instance Ord MasterRecord where
    compare = comparing mchrom
           <> comparing mtargetstart
           <> comparing mtargetend

defaultMasterRec = MasterRecord { mchrom = NONE
                                , mtargetstart = 0
                                , mtargetend = 0
                                , mampname = "EMPTY"
                                , mforstart = 0
                                , mforend = 0
                                , mrevstart = 0
                                , mrevend = 0
                                , mforname = "EMPTY"
                                , mrevname = "EMPTY"
                                , mforseq = "EMPTY"
                                , mrevseq = "EMPTY"
                                }

-- 180226 record type to hold run stats
data RunStats = RunStats { alnsTotal :: Integer
                         , alnsMapped :: Integer
                         , alnsTrimd :: Integer
                         , alnsTrimdToZero :: Integer
                         , trimmedPct :: Double
                         , mappedPct :: Double
                         } deriving (Eq, Show, Read)

-- 170919 add MidFamily record type to group [AlignedRead] by mID
data MidFamily = MidFamily { chrom :: UChr
                           , alnstart :: Integer
                           , alnend :: Integer
                           , midseq :: B.ByteString
                           , alns :: [AlignedRead]
                           , alncount :: Integer
                           , txposedcount :: Integer
                           , postxposdcnt :: Integer
                           , negtxposdcnt :: Integer
                           , targetbed :: BedRecord
                           } deriving (Show, Eq)

instance Ord MidFamily where
    compare = comparing chrom
           <> comparing alnstart
           <> comparing alnend
           <> comparing alncount

defaultMidFam = MidFamily { chrom = NONE
                          , alnstart = 0
                          , alnend = 0
                          , midseq = "EMPTY"
                          , alns = []
                          , alncount = 0
                          , txposedcount = 0
                          , postxposdcnt = 0
                          , negtxposdcnt = 0
                          , targetbed = defaultBed
                          }

type Header = [B.ByteString]
type Alignments = [AlignedRead]
type BED = V.Vector BedRecord
type SAM = (Header, Alignments)
type CMap = M.Map UChr BedMap
type BedMap = I.IntMap BedRecord
type CigarMap = [(Integer, B.ByteString)]

-- Internal sum type for all chromosome names in both UCSC and ensembl naming
-- schemes.
-- NOTE: only GRCh37 non-canonical chromosomes recognized this version
data UChr = C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 | C12 | C13
          | C14 | C15 | C16 | C17 | C18 | C19 | C20 | C21 | C22 | CX | CY | CMT
          | Chr1 | Chr2 | Chr3 | Chr4 | Chr5 | Chr6 | Chr7 | Chr8 | Chr9
          | Chr10 | Chr11 | Chr12 | Chr13 | Chr14 | Chr15 | Chr16 | Chr17
          | Chr18 | Chr19 | Chr20 | Chr21 | Chr22 | ChrX | ChrY | ChrM
          | GL000207P1
          | GL000226P1
          | GL000229P1
          | GL000231P1
          | GL000210P1
          | GL000239P1
          | GL000235P1
          | GL000201P1
          | GL000247P1
          | GL000245P1
          | GL000197P1
          | GL000203P1
          | GL000246P1
          | GL000249P1
          | GL000196P1
          | GL000248P1
          | GL000244P1
          | GL000238P1
          | GL000202P1
          | GL000234P1
          | GL000232P1
          | GL000206P1
          | GL000240P1
          | GL000236P1
          | GL000241P1
          | GL000243P1
          | GL000242P1
          | GL000230P1
          | GL000237P1
          | GL000233P1
          | GL000204P1
          | GL000198P1
          | GL000208P1
          | GL000191P1
          | GL000227P1
          | GL000228P1
          | GL000214P1
          | GL000221P1
          | GL000209P1
          | GL000218P1
          | GL000220P1
          | GL000213P1
          | GL000211P1
          | GL000199P1
          | GL000217P1
          | GL000216P1
          | GL000215P1
          | GL000205P1
          | GL000219P1
          | GL000224P1
          | GL000223P1
          | GL000195P1
          | GL000212P1
          | GL000222P1
          | GL000200P1
          | GL000193P1
          | GL000194P1
          | GL000225P1
          | GL000192P1
          | NC_007605
          | NONE
            deriving (Eq, Ord, Generic)

-- 170508 type for tracking chromosome name type
data ChromNameFormat = GRC | UCSC
                        deriving (Show, Eq, Ord, Bounded, Enum, Read)

instance Show UChr where
    show C1 = "1"
    show C2 = "2"
    show C3 = "3"
    show C4 = "4"
    show C5 = "5"
    show C6 = "6"
    show C7 = "7"
    show C8 = "8"
    show C9 = "9"
    show C10 = "10"
    show C11 = "11"
    show C12 = "12"
    show C13 = "13"
    show C14 = "14"
    show C15 = "15"
    show C16 = "16"
    show C17 = "17"
    show C18 = "18"
    show C19 = "19"
    show C20 = "20"
    show C21 = "21"
    show C22 = "22"
    show CX = "X"
    show CY = "Y"
    show CMT = "MT"
    show GL000207P1 = "GL000207.1"
    show GL000226P1 = "GL000226.1"
    show GL000229P1 = "GL000229.1"
    show GL000231P1 = "GL000231.1"
    show GL000210P1 = "GL000210.1"
    show GL000239P1 = "GL000239.1"
    show GL000235P1 = "GL000235.1"
    show GL000201P1 = "GL000201.1"
    show GL000247P1 = "GL000247.1"
    show GL000245P1 = "GL000245.1"
    show GL000197P1 = "GL000197.1"
    show GL000203P1 = "GL000203.1"
    show GL000246P1 = "GL000246.1"
    show GL000249P1 = "GL000249.1"
    show GL000196P1 = "GL000196.1"
    show GL000248P1 = "GL000248.1"
    show GL000244P1 = "GL000244.1"
    show GL000238P1 = "GL000238.1"
    show GL000202P1 = "GL000202.1"
    show GL000234P1 = "GL000234.1"
    show GL000232P1 = "GL000232.1"
    show GL000206P1 = "GL000206.1"
    show GL000240P1 = "GL000240.1"
    show GL000236P1 = "GL000236.1"
    show GL000241P1 = "GL000241.1"
    show GL000243P1 = "GL000243.1"
    show GL000242P1 = "GL000242.1"
    show GL000230P1 = "GL000230.1"
    show GL000237P1 = "GL000237.1"
    show GL000233P1 = "GL000233.1"
    show GL000204P1 = "GL000204.1"
    show GL000198P1 = "GL000198.1"
    show GL000208P1 = "GL000208.1"
    show GL000191P1 = "GL000191.1"
    show GL000227P1 = "GL000227.1"
    show GL000228P1 = "GL000228.1"
    show GL000214P1 = "GL000214.1"
    show GL000221P1 = "GL000221.1"
    show GL000209P1 = "GL000209.1"
    show GL000218P1 = "GL000218.1"
    show GL000220P1 = "GL000220.1"
    show GL000213P1 = "GL000213.1"
    show GL000211P1 = "GL000211.1"
    show GL000199P1 = "GL000199.1"
    show GL000217P1 = "GL000217.1"
    show GL000216P1 = "GL000216.1"
    show GL000215P1 = "GL000215.1"
    show GL000205P1 = "GL000205.1"
    show GL000219P1 = "GL000219.1"
    show GL000224P1 = "GL000224.1"
    show GL000223P1 = "GL000223.1"
    show GL000195P1 = "GL000195.1"
    show GL000212P1 = "GL000212.1"
    show GL000222P1 = "GL000222.1"
    show GL000200P1 = "GL000200.1"
    show GL000193P1 = "GL000193.1"
    show GL000194P1 = "GL000194.1"
    show GL000225P1 = "GL000225.1"
    show GL000192P1 = "GL000192.1"
    show NC_007605 = "NC_007605"
    show Chr1 = "chr1"
    show Chr2 = "chr2"
    show Chr3 = "chr3"
    show Chr4 = "chr4"
    show Chr5 = "chr5"
    show Chr6 = "chr6"
    show Chr7 = "chr7"
    show Chr8 = "chr8"
    show Chr9 = "chr9"
    show Chr10 = "chr10"
    show Chr11 = "chr11"
    show Chr12 = "chr12"
    show Chr13 = "chr13"
    show Chr14 = "chr14"
    show Chr15 = "chr15"
    show Chr16 = "chr16"
    show Chr17 = "chr17"
    show Chr18 = "chr18"
    show Chr19 = "chr19"
    show Chr20 = "chr20"
    show Chr21 = "chr21"
    show Chr22 = "chr22"
    show ChrX = "chrX"
    show ChrY = "chrY"
    show ChrM = "chrM"
    show NONE = "*"

{--
-- 170508 add alternate "show" function to print chromosome names based on
-- config flag setting
showChrom :: ChromNameFormat -> UChr -> String
showChrom fmt chr = case (fmt, chr) of
    (GRC,Chr1)   -> "1"
    (GRC,Chr2)   -> "2"
    (GRC,Chr3)   -> "3"
    (GRC,Chr4)   -> "4"
    (GRC,Chr5)   -> "5"
    (GRC,Chr6)   -> "6"
    (GRC,Chr7)   -> "7"
    (GRC,Chr8)   -> "8"
    (GRC,Chr9)   -> "9"
    (GRC,Chr10)  -> "10"
    (GRC,Chr11)  -> "11"
    (GRC,Chr12)  -> "12"
    (GRC,Chr13)  -> "13"
    (GRC,Chr14)  -> "14"
    (GRC,Chr15)  -> "15"
    (GRC,Chr16)  -> "16"
    (GRC,Chr17)  -> "17"
    (GRC,Chr18)  -> "18"
    (GRC,Chr19)  -> "19"
    (GRC,Chr20)  -> "20"
    (GRC,Chr21)  -> "21"
    (GRC,Chr22)  -> "22"
    (GRC,ChrX)   -> "X"
    (GRC,ChrY)   -> "Y"
    (GRC,ChrM)   -> "MT"
    (UCSC,Chr1)  -> "chr1"
    (UCSC,Chr2)  -> "chr2"
    (UCSC,Chr3)  -> "chr3"
    (UCSC,Chr4)  -> "chr4"
    (UCSC,Chr5)  -> "chr5"
    (UCSC,Chr6)  -> "chr6"
    (UCSC,Chr7)  -> "chr7"
    (UCSC,Chr8)  -> "chr8"
    (UCSC,Chr9)  -> "chr9"
    (UCSC,Chr10) -> "chr10"
    (UCSC,Chr11) -> "chr11"
    (UCSC,Chr12) -> "chr12"
    (UCSC,Chr13) -> "chr13"
    (UCSC,Chr14) -> "chr14"
    (UCSC,Chr15) -> "chr15"
    (UCSC,Chr16) -> "chr16"
    (UCSC,Chr17) -> "chr17"
    (UCSC,Chr18) -> "chr18"
    (UCSC,Chr19) -> "chr19"
    (UCSC,Chr20) -> "chr20"
    (UCSC,Chr21) -> "chr21"
    (UCSC,Chr22) -> "chr22"
    (UCSC,ChrX)  -> "chrX"
    (UCSC,ChrY)  -> "chrY"
    (UCSC,ChrM)  -> "chrM"
    otherwise    -> "*"
--}

-- 170510 try parsing either GRC or UCSC chromosome names
uchrparser :: A.Parser UChr
uchrparser =      (A.string "chr10" >> return Chr10)
        <|>  (A.string "chr11" >> return Chr11)
        <|>  (A.string "chr12" >> return Chr12)
        <|>  (A.string "chr13" >> return Chr13)
        <|>  (A.string "chr14" >> return Chr14)
        <|>  (A.string "chr15" >> return Chr15)
        <|>  (A.string "chr16" >> return Chr16)
        <|>  (A.string "chr17" >> return Chr17)
        <|>  (A.string "chr18" >> return Chr18)
        <|>  (A.string "chr19" >> return Chr19)
        <|>  (A.string "chr20" >> return Chr20)
        <|>  (A.string "chr21" >> return Chr21)
        <|>  (A.string "chr22" >> return Chr22)
        <|>  (A.string "chr1" >> return Chr1)
        <|>  (A.string "chr2" >> return Chr2)
        <|>  (A.string "chr3" >> return Chr3)
        <|>  (A.string "chr4" >> return Chr4)
        <|>  (A.string "chr5" >> return Chr5)
        <|>  (A.string "chr6" >> return Chr6)
        <|>  (A.string "chr7" >> return Chr7)
        <|>  (A.string "chr8" >> return Chr8)
        <|>  (A.string "chr9" >> return Chr9)
        <|>  (A.string "chrX" >> return ChrX)
        <|>  (A.string "chrY" >> return ChrY)
        <|>  (A.string "chrM" >> return ChrM)
        <|>  (A.string "10" >> return C10)
        <|>  (A.string "11" >> return C11)
        <|>  (A.string "12" >> return C12)
        <|>  (A.string "13" >> return C13)
        <|>  (A.string "14" >> return C14)
        <|>  (A.string "15" >> return C15)
        <|>  (A.string "16" >> return C16)
        <|>  (A.string "17" >> return C17)
        <|>  (A.string "18" >> return C18)
        <|>  (A.string "19" >> return C19)
        <|>  (A.string "20" >> return C20)
        <|>  (A.string "21" >> return C21)
        <|>  (A.string "22" >> return C22)
        <|>  (A.string "1" >> return C1)
        <|>  (A.string "2" >> return C2)
        <|>  (A.string "3" >> return C3)
        <|>  (A.string "4" >> return C4)
        <|>  (A.string "5" >> return C5)
        <|>  (A.string "6" >> return C6)
        <|>  (A.string "7" >> return C7)
        <|>  (A.string "8" >> return C8)
        <|>  (A.string "9" >> return C9)
        <|>  (A.string "X" >> return CX)
        <|>  (A.string "Y" >> return CY)
        <|>  (A.string "MT" >> return CMT)
        <|>  (A.string "GL000207.1" >> return GL000207P1)
        <|>  (A.string "GL000226.1" >> return GL000226P1)
        <|>  (A.string "GL000229.1" >> return GL000229P1)
        <|>  (A.string "GL000231.1" >> return GL000231P1)
        <|>  (A.string "GL000210.1" >> return GL000210P1)
        <|>  (A.string "GL000239.1" >> return GL000239P1)
        <|>  (A.string "GL000235.1" >> return GL000235P1)
        <|>  (A.string "GL000201.1" >> return GL000201P1)
        <|>  (A.string "GL000247.1" >> return GL000247P1)
        <|>  (A.string "GL000245.1" >> return GL000245P1)
        <|>  (A.string "GL000197.1" >> return GL000197P1)
        <|>  (A.string "GL000203.1" >> return GL000203P1)
        <|>  (A.string "GL000246.1" >> return GL000246P1)
        <|>  (A.string "GL000249.1" >> return GL000249P1)
        <|>  (A.string "GL000196.1" >> return GL000196P1)
        <|>  (A.string "GL000248.1" >> return GL000248P1)
        <|>  (A.string "GL000244.1" >> return GL000244P1)
        <|>  (A.string "GL000238.1" >> return GL000238P1)
        <|>  (A.string "GL000202.1" >> return GL000202P1)
        <|>  (A.string "GL000234.1" >> return GL000234P1)
        <|>  (A.string "GL000232.1" >> return GL000232P1)
        <|>  (A.string "GL000206.1" >> return GL000206P1)
        <|>  (A.string "GL000240.1" >> return GL000240P1)
        <|>  (A.string "GL000236.1" >> return GL000236P1)
        <|>  (A.string "GL000241.1" >> return GL000241P1)
        <|>  (A.string "GL000243.1" >> return GL000243P1)
        <|>  (A.string "GL000242.1" >> return GL000242P1)
        <|>  (A.string "GL000230.1" >> return GL000230P1)
        <|>  (A.string "GL000237.1" >> return GL000237P1)
        <|>  (A.string "GL000233.1" >> return GL000233P1)
        <|>  (A.string "GL000204.1" >> return GL000204P1)
        <|>  (A.string "GL000198.1" >> return GL000198P1)
        <|>  (A.string "GL000208.1" >> return GL000208P1)
        <|>  (A.string "GL000191.1" >> return GL000191P1)
        <|>  (A.string "GL000227.1" >> return GL000227P1)
        <|>  (A.string "GL000228.1" >> return GL000228P1)
        <|>  (A.string "GL000214.1" >> return GL000214P1)
        <|>  (A.string "GL000221.1" >> return GL000221P1)
        <|>  (A.string "GL000209.1" >> return GL000209P1)
        <|>  (A.string "GL000218.1" >> return GL000218P1)
        <|>  (A.string "GL000220.1" >> return GL000220P1)
        <|>  (A.string "GL000213.1" >> return GL000213P1)
        <|>  (A.string "GL000211.1" >> return GL000211P1)
        <|>  (A.string "GL000199.1" >> return GL000199P1)
        <|>  (A.string "GL000217.1" >> return GL000217P1)
        <|>  (A.string "GL000216.1" >> return GL000216P1)
        <|>  (A.string "GL000215.1" >> return GL000215P1)
        <|>  (A.string "GL000205.1" >> return GL000205P1)
        <|>  (A.string "GL000219.1" >> return GL000219P1)
        <|>  (A.string "GL000224.1" >> return GL000224P1)
        <|>  (A.string "GL000223.1" >> return GL000223P1)
        <|>  (A.string "GL000195.1" >> return GL000195P1)
        <|>  (A.string "GL000212.1" >> return GL000212P1)
        <|>  (A.string "GL000222.1" >> return GL000222P1)
        <|>  (A.string "GL000200.1" >> return GL000200P1)
        <|>  (A.string "GL000193.1" >> return GL000193P1)
        <|>  (A.string "GL000194.1" >> return GL000194P1)
        <|>  (A.string "GL000225.1" >> return GL000225P1)
        <|>  (A.string "GL000192.1" >> return GL000192P1)
        <|>  (A.string "NC_007605" >> return NC_007605)
        <|>  (A.string "*" >> return NONE)
        -- <|>  otherchromp -- 181024

{--
otherchromp :: A.Parser UChr
otherchromp = do
    c <- A.takeTill A.isSpace
    return $ Other c
--}

-- 180206 include option for providing primer coords in BED or BEDPE format
-- in addition to master file
optargs :: Parser Opts
optargs = Opts
    <$> switch
        ( short 'b'
       <> long "bedpe"
       <> help "add this switch to use BEDPE coordinate input format (default format is master file)"
        )
    <*> argument str (metavar "PRIMER_COORDS_INFILE")
    <*> argument str (metavar "SAM_INFILE")
    <*> argument str (metavar "OUTPUT_SAM_FILENAME")

-- record to store command line arguments
data Opts = Opts { bedpeformat :: Bool
                 , incoordsfile :: String
                 , insamfile :: String
                 , outfilename :: String
                 } deriving (Show, Eq)

-- 170927 parse master file for primer and target intervals
masterparser :: A.Parser MasterRecord
masterparser = do
    tbed <- bedrecparser
    A.skipSpace
    fpstart <- A.decimal
    A.skipSpace
    fpend <- A.decimal
    A.skipSpace
    fpname <- txtfieldp
    A.skipSpace
    rpstart <- A.decimal
    A.skipSpace
    rpend <- A.decimal
    A.skipSpace
    rpname <- txtfieldp
    A.skipSpace
    fseq <- seqp
    A.skipSpace
    rseq <- seqp
    return  MasterRecord { mchrom = bedchr tbed
                         , mtargetstart = bedstart tbed
                         , mtargetend = bedend tbed
                         , mampname = bedname tbed
                         , mforstart = fpstart
                         , mforend = fpend
                         , mrevstart = rpstart
                         , mrevend = rpend
                         , mforname = fpname
                         , mrevname = rpname
                         , mforseq = fseq
                         , mrevseq = rseq
                         }

bedrecparser :: A.Parser BedRecord
bedrecparser = do
    chr <- uchrparser
    A.skipSpace
    bstart <- A.decimal
    A.skipSpace
    bend <- A.decimal
    A.skipSpace
    bname <- txtfieldp
    let bedrec = defaultBed { bedchr = chr
                            , bedstart = bstart
                            , bedend = bend
                            , bedname = bname
                            }
    return bedrec

-- 20180206
bedPEparser :: A.Parser BEDPE
bedPEparser = do
    c1 <- uchrparser
    A.skipSpace
    s1 <- A.decimal
    A.skipSpace
    e1 <- A.decimal
    A.skipSpace
    c2 <- uchrparser
    A.skipSpace
    s2 <- A.decimal
    A.skipSpace
    e2 <- A.decimal
    bname <- txtfieldp
    -- ignore any additional columns present
    return $ BEDPE c1 s1 e1 c2 s2 e2 bname

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- read SAM file and parse header and alignments, returning both in tuple
-- NOTE: as of October 2017 these functions are deprecated due to use of the
-- Conduit library for lower memory usage. They remain in the library in the
-- event they are needed for future use.
------------------------------------------------------------------------------
------------------------------------------------------------------------------

-- 180321 parse SAM file into PairedAln records to group alignments by read name
-- NOTE: first approach will be to filter for alignments intersecting one or
-- more primer intervals, and only running those alignments through the trimming.
readSAMtoPairedAlns :: FilePath -> IO (AlignedRead, [PairedAln])
readSAMtoPairedAlns fp = do
    parsedAlns <- readSAMnameset fp
    let hdr = safegetheader parsedAlns
        pairsets = alnsToPairedAln <$> (tail parsedAlns)
    return (hdr, pairsets)

alnsToPairedAln :: [AlignedRead] -> PairedAln
alnsToPairedAln [] = defaultPairedAln
alnsToPairedAln as =
    let (r1s, r2s) = partition read1filter as
        (r1pl, r1secs) = partition primaryR1filter r1s
        (r2pl, r2secs) = partition primaryR2filter r2s
        r1p = headsafeAln r1pl
        r2p = headsafeAln r2pl
    in PairedAln r1p r2p r1secs r2secs

-- 180321 non-conduit read name-set SAM file reading function
readSAMnameset :: FilePath -> IO [[AlignedRead]]
readSAMnameset fp = do
    txt <- B.readFile fp
    return $ parseReadsetsFromSAM txt

-- 180212 non-conduit latest SAM parsing function for debugging
readSAM :: FilePath -> IO [AlignedRead]
readSAM fp = do
    bslines <- B.lines <$> B.readFile fp
    let parsedAlns = (A.parseOnly (hdralnparser <|> alnparser)) <$> bslines
    return $ U.rights parsedAlns

-- NOTE: changed to make use of U.rights to leave out parsing failures
parseAlns :: [B.ByteString] -> Alignments
parseAlns as = U.rights $ (A.parseOnly alnparser <$> as)

parseAln as = A.parseOnly alnparser as

-- 180329 parse to PairedAln directly
parsePairedAlnsFromSAM :: B.ByteString -> Either String [PairedAln]
parsePairedAlnsFromSAM bs = A.parseOnly parsePairedAlns bs

parsePairedAlnsOrHdr = A.many1 (hdralnparserEOL <|> pairedalnparser)

parsePairedAlns = do
    h <- hdralnparserEOL
    pairsets <- A.many1 pairedalnparser
    return $ h : pairsets

pairedalnparser = do
    r1 <- alnparser
    A.endOfLine
    let setname = qname r1
    secrs <- A.many1 ((secalnp setname) <* A.endOfLine)
    return $ alnsToPairedAln (r1 : secrs)

-- 180321 parse SAM input as continuous text (do not split into lines) so that
-- the lists of reads (grouped by read name) can be organized into paired alignments
parseReadsetsFromSAM :: B.ByteString -> [[AlignedRead]]
parseReadsetsFromSAM bs = U.fromRight [] $ A.parseOnly parsePairedAlns' bs

parsePairedAlns' = do
    h <- hdralnparserEOL'
    rsets <- A.many1 alignmentsetparser
    return $ [h] : rsets

-- 180321 try out a sort of backref for unknown number of secondary alignments
alignmentsetparser = do
    r1 <- alnparser
    A.endOfLine
    let setname = qname r1
    secrs <- A.many1 ((secalnp setname) <* A.endOfLine)
    return $ r1 : secrs

secalnp :: B.ByteString -> A.Parser AlignedRead
secalnp setqname = do
    qn <- A.string setqname
    A.space
    f <- A.decimal
    A.space
    chr <- uchrparser
    A.space
    p <- A.decimal
    A.space
    mpscore <- A.decimal
    A.space
    cig <- txtfieldp
    A.space
    nchr <- txtfieldp
    A.space
    pn <- A.decimal
    A.space
    tl <- A.double
    A.space
    seq <- txtfieldp
    A.space
    qual <- txtfieldp
    A.space
    optfs <- optfieldstotalp -- optfieldsp
    -- A.endOfLine -- 180321 for parsing alignment pairs
    let flag = f
        strand = case testBit flag 4 of
            True  -> B.pack "-"
            False -> B.pack "+"
        prd = testBit flag 0
        mapd = not $ testBit flag 2
        tln = floor $ tl
        cigm = exrights $ parseCigar cig
        end = (sumMatches cigm) + p -- 180223 NOTE: previous CIGAR accounting inverted "D" and "I" (!)
        -- optfieldstr = B.intercalate "\t" optfs -- B.pack optfs
        midstr = parsemIDstring optfs -- optfieldstr
        a = defaultAlignment { qname = qn
                             , flag = f
                             , rname = chr
                             , pos = p - 1 --  SAM to BED/BAM numbering
                             , trimdpos = p - 1 --  SAM to BED/BAM numbering
                             , endpos = end - 1 --  SAM to BED/BAM numbering
                             , trimdendpos = end - 1 --  SAM to BED/BAM numbering
                             , mapqual = mpscore
                             , cigar = cig
                             , cigmap = cigm
                             , rnext = nchr
                             , pnext = pn - 1 -- SAM to BED/BAM
                             , tlen = tln
                             , refseq = seq
                             , basequal = qual
                             , optfields = optfs -- optfieldstr
                             , strand = strand
                             , paired = prd
                             , mapped = mapd
                             , mid = midstr
                             }
    return a


readBEDPE :: FilePath -> IO [BEDPE]
readBEDPE fp = do
    flines <- B.lines <$> B.readFile fp
    let nr = length flines
        bs = parseBEDPE <$> flines -- [BEDPE]
        lks = [0..] :: [Int]
        bmap = M.fromList $ zip lks bs
        (succm, failm) = M.partition U.isRight bmap
        succs = U.rights $ snd <$> M.toAscList succm
        faildlineixs = M.keys failm
        bcnt = length succs
        failcnt = length faildlineixs
        parsestatus = parsechkBED nr bcnt faildlineixs
        parsefaildLines = B.unlines $ (flines !!) <$> faildlineixs
    putStrLn parsestatus
    B.writeFile "bedPEparsefails.log" parsefaildLines
    writeFile "primer_BEDPE_parsing.log" parsestatus
    return succs

parseBED bs = A.parseOnly bedrecparser bs

parseBEDPE bs = A.parseOnly bedPEparser bs

-- 170510 write parse fail lines to log file (don't print line nums to stdout)
parsechkSAM :: Int -> Int -> [Int] -> String
parsechkSAM numrec numsucc failedlines
    | numrec == numsucc = reportOK
    | otherwise = reportFails
        where reportOK = "all " ++ (show numrec)
                                ++ " SAM alignments parsed successfully.\n"
              parsediff = length failedlines
              --flinestrs = show <$> failedlines
              --fstring = intercalate "\n" flinestrs
              reportFails = "WARNING: "
                             ++ (show parsediff)
                             ++ " of " ++ (show numrec)
                             ++ " SAM alignments failed to parse (see samparsefails.log)"
                             ++ "\n"
                             -- ++ fstring ++ "\n"

parsechkBED :: Int -> Int -> [Int] -> String
parsechkBED numrec numsucc failedlines
    | numrec == numsucc = reportOK
    | otherwise = reportFails
        where reportOK = "all " ++ (show numrec)
                                ++ " primer BEDPE records parsed successfully.\n"
              parsediff = length failedlines
              -- flinestrs = show <$> failedlines
              -- fstring = intercalate "\n" flinestrs
              reportFails = "WARNING: "
                             ++ (show parsediff)
                             ++ " of " ++ (show numrec)
                             ++ " primer BEDPE records failed to parse (see bedparsefails.log)"
                             ++ "\n"


-- 170926 read and parse master file
getMasterFile :: FilePath -> IO [MasterRecord]
getMasterFile fp = do
    mlines <- B.lines <$> B.readFile fp
    let nr = length mlines
        mrecs = (A.parseOnly masterparser) <$> mlines
        lks = [0..] :: [Int]
        mrecmap = M.fromList $ zip lks mrecs
        (succm, failm) = M.partition U.isRight mrecmap
        succs = U.rights $ snd <$> M.toAscList succm
        failedlinenums = M.keys failm
        mcnt = length succs
        failcnt = length failedlinenums
        parseStatus = parsechkMaster nr mcnt failedlinenums
        parsefailrecs = B.unlines $ (mlines !!) <$> failedlinenums
    putStrLn parseStatus
    B.writeFile "masterparsefails.log" $ parsefailrecs
    return succs

parsechkMaster :: Int -> Int -> [Int] -> String
parsechkMaster numrec numsucc failedlines
    | numrec == numsucc = reportOK
    | otherwise = reportFails
        where reportOK = "all " ++ (show numrec)
                                ++ " master records parsed successfully.\n"
              parsediff = length failedlines
              reportFails = "WARNING: "
                             ++ (show parsediff)
                             ++ " of " ++ (show numrec)
                             ++ " master records failed to parse (see masterparsefails.log)"
                             ++ "\n"

masterToFPrimerBED :: [MasterRecord] -> BED
masterToFPrimerBED mrecs = V.fromList
                         $ getFPrimerBEDfromMaster <$> (sort mrecs)

masterToRPrimerBED :: [MasterRecord] -> BED
masterToRPrimerBED mrecs = V.fromList
                         $ getRPrimerBEDfromMaster <$> (sort mrecs)

masterToTargetBED :: [MasterRecord] -> BED
masterToTargetBED mrecs = V.fromList
                         $ getTargetBEDfromMaster <$> (sort mrecs)

masterRecToTargetBED :: MasterRecord -> BedRecord
masterRecToTargetBED m = defaultBed { bedchr = mchrom m
                                    , bedstart = mtargetstart m
                                    , bedend = mtargetend m
                                    , bedname = mampname m
                                    }

getFPrimerBEDfromMaster :: MasterRecord -> BedRecord
getFPrimerBEDfromMaster mr =
    let c = mchrom mr
        fpstart = mforstart mr
        fpend = mforend mr
        fpname = mforname mr
    in defaultBed { bedchr = c
                  , bedstart = fpstart
                  , bedend = fpend
                  , bedname = fpname
                  }

getRPrimerBEDfromMaster :: MasterRecord -> BedRecord
getRPrimerBEDfromMaster mr =
    let c = mchrom mr
        rpstart = mrevstart mr
        rpend = mrevend mr
        rpname = mrevname mr
    in defaultBed { bedchr = c
                  , bedstart = rpstart
                  , bedend = rpend
                  , bedname = rpname
                  }

getTargetBEDfromMaster :: MasterRecord -> BedRecord
getTargetBEDfromMaster mr =
    let c = mchrom mr
        tstart = mtargetstart mr
        tend = mtargetend mr
        tname = mampname mr
    in defaultBed { bedchr = c
                  , bedstart = tstart
                  , bedend = tend
                  , bedname = tname
                  }

-- new alignment parser which allows parse failures to be dropped
-- from Alignments
alnparser :: A.Parser AlignedRead
alnparser = do
    qn <- txtfieldp
    A.space
    f <- A.decimal
    A.space
    chr <- uchrparser
    A.space
    p <- A.decimal
    A.space
    mpscore <- A.decimal
    A.space
    cig <- txtfieldp
    A.space
    nchr <- txtfieldp
    A.space
    pn <- A.decimal
    A.space
    tl <- A.double
    A.space
    seq <- txtfieldp
    A.space
    qual <- txtfieldp
    A.space
    optfs <- optfieldstotalp -- optfieldsp
    -- A.endOfLine -- 180321 for parsing alignment pairs
    let flag = f
        strand = case testBit flag 4 of
            True  -> B.pack "-"
            False -> B.pack "+"
        prd = testBit flag 0
        mapd = not $ testBit flag 2
        tln = floor $ tl
        cigm = exrights $ parseCigar cig
        end = (sumMatches cigm) + p -- 180223 NOTE: previous CIGAR accounting inverted "D" and "I" (!)
        -- optfieldstr = B.intercalate "\t" optfs -- B.pack optfs
        midstr = parsemIDstring optfs -- optfieldstr
        a = defaultAlignment { qname = qn
                             , flag = f
                             , rname = chr
                             , pos = p - 1 --  SAM to BED/BAM numbering
                             , trimdpos = p - 1 --  SAM to BED/BAM numbering
                             , endpos = end - 1 --  SAM to BED/BAM numbering
                             , trimdendpos = end - 1 --  SAM to BED/BAM numbering
                             , mapqual = mpscore
                             , cigar = cig
                             , cigmap = cigm
                             , rnext = nchr
                             , pnext = pn - 1 -- SAM to BED/BAM
                             , tlen = tln
                             , refseq = seq
                             , basequal = qual
                             , optfields = optfs -- optfieldstr
                             , strand = strand
                             , paired = prd
                             , mapped = mapd
                             , mid = midstr
                             }
    return a

-- 180329 keep homogenous type throughout conduit stream
hdralnparserEOL :: A.Parser PairedAln
hdralnparserEOL = do
    hlines <- A.many1 samhdrparserEOL
    let hdraln = defaultAlignment { headerstrings = hlines
                                  , isheader = True
                                  , qname = "HEADERLINE"
                                  }
    return $ defaultPairedAln { r1prim = hdraln }

hdralnparserEOL' :: A.Parser AlignedRead
hdralnparserEOL' = do
    hlines <- A.many1 samhdrparserEOL
    return $ defaultAlignment { headerstrings = hlines
                              , isheader = True
                              , qname = "HEADERLINE"
                              }

-- 171017 header parser more suited for conduit and parsing header into AlignedRead
hdralnparser :: A.Parser AlignedRead
hdralnparser = do
    hlines <- A.many1 samhdrparser
    return $ defaultAlignment { headerstrings = hlines
                              , isheader = True
                              , qname = "HEADERLINE"
                              }

-- SAM header chromosome names parser (parses chromosome names ONLY for now)
-- NOTE: for simplicity, only canonical chromosome names are parsed at this time.
hdrchromp :: A.Parser UChr
hdrchromp = do
    A.string "@SQ"
    A.skipSpace
    A.string "SN:"
    c <- uchrparser
    A.many1 A.anyChar
    return c

samhdrparser :: A.Parser B.ByteString
samhdrparser = do
    atprefix <- A.char '@'
    restofhdrline <- A.takeTill (A.inClass "\n\r")
    let hdrln = B.append "@" restofhdrline
    return hdrln

-- 171204 parser for header line with newline char present
samhdrparserEOL :: A.Parser B.ByteString
samhdrparserEOL = do
    atprefix <- A.char '@'
    fields <- A.takeTill (A.inClass "\n\r")
    A.endOfLine
    let hdrln = B.append "@" fields
    return hdrln

-- 180409 generic filter function over PairedAln
anyPairedAln :: (AlignedRead -> Bool) -> PairedAln -> Bool
anyPairedAln b p
    | contains = True
    | otherwise = False
        where contains = r1pred || r2pred || r1secspred || r2secspred
              r1pred = (b $ r1prim p)
              r2pred = (b $ r2prim p)
              r1secspred = (any b $ r1secs p)
              r2secspred = (any b $ r2secs p)

allPairedAln :: (AlignedRead -> Bool) -> PairedAln -> Bool
allPairedAln b p
    | contains = True
    | otherwise = False
        where contains = r1pred && r2pred && r1secspred && r2secspred
              r1pred = (b $ r1prim p)
              r2pred = (b $ r2prim p)
              r1secspred = (any b $ r1secs p)
              r2secspred = (any b $ r2secs p)

-- 180321 safely get header AlignedRead from [[AlignedRead]]
safegetheader :: [[AlignedRead]] -> AlignedRead
safegetheader as
    | null $ join as = error "[!] SAM header parse error: check header in input SAM file" -- defaultAlignment
    | isheader hdr   = hdr
    | otherwise      = error "[!] SAM header missing: check header in input SAM file" -- defaultAlignment
        where hdr = head $ head as

headsafeAln :: [AlignedRead] -> AlignedRead
headsafeAln as
    | null as = error "[!] SAM header missing: check header in input SAM file"
    | otherwise = head as

-- 171204 Attoparsec version of getRight (below)
-- use of Conduit to read input alignment file precludes use of "rights" in Data.Either
-- unless/until we come up with a more idiomatic implementation of the parsing step.
rightOrDefault :: Either String AlignedRead -> AlignedRead
rightOrDefault e = case e of
    Left _ -> defaultAlignment
    Right a -> a

rightOrDefaultPaird :: Either CA.ParseError (CA.PositionRange, [PairedAln])
                    -> [PairedAln]
rightOrDefaultPaird e = case e of
    Left _ -> [] -- defaultPairedAln
    Right a -> snd a

-- 171017 extract parse success from Either (as singleton, for conduit)
getRight :: Either t (a, AlignedRead) -> AlignedRead
getRight e = case e of
    Left _  -> defaultAlignment
    Right x -> snd x

-- print parse failures to log file (try reporting attoparsec output to start)
-- getRightM :: (P.MonadIO m, Show t) => Either t (a, AlignedRead) -> IO AlignedRead
getRightM ev = case ev of
    Left e  -> do
        -- let epos = errorPosition e
        liftIO $ writeFile "samparsefails.log" $ show e
        return defaultAlignment
    Right a -> return $ snd a

-- cigar parsers
starcigarP = A.string "*" >> return [(0, star)]
    where star = "*" :: B.ByteString

cigarP = A.many1 ((,) <$> A.decimal <*> (A.takeWhile (A.inClass "MIDNSHP=X")))

parseCigar = A.parseOnly $ cigarP <|> starcigarP

-- used on single cigar ByteString
getAlignedLength cigcol = sumMatches ciglist
    where ciglist = head $ U.rights $ cig
          cig = (parseCigar cigcol) : []

sumMatches :: [(Integer, B.ByteString)] -> Integer
sumMatches cigs = sum [ x | (x, y) <- cigs, y == "M" || y == "D" ]

sumSeqMatches :: [(Integer, B.ByteString)] -> Integer
sumSeqMatches cigs = sum [ x | (x, y) <- cigs, y == "M"
                                            || y == "I"
                                            || y == "S" ]

sumRefMatches :: [(Integer, B.ByteString)] -> Integer
sumRefMatches cigs = sum [ x | (x, y) <- cigs, y == "M" || y == "D" ]

sumSoftClipCigOps :: [(Integer, B.ByteString)] -> Int
sumSoftClipCigOps cigs = genericLength [ o | (_, o) <- cigs, o == "M"
                                                          || o == "I" ]

-- 180212
nomapCigToNomapRname :: B.ByteString -> UChr -> UChr
nomapCigToNomapRname c rname
    | c == "*"  = NONE
    | otherwise = rname

-- 180223 check if cigar "length" equal to alignment length (refseq)


-- 180223 check that trimmed cigar string still matches read length
checkcigseqlen :: AlignedRead -> Bool
checkcigseqlen a
    | (trimdcigar a == "*") || (cigar a == "*") = True
    | cigmatchlen == refseqlen = True
    | otherwise = False
        where cigmatchlen = sumSeqMatches tcmap
              tcmap = mapcig $ trimdcigar a
              refseqlen = genericLength $ B.unpack $ refseq a

getcigseqdiff :: AlignedRead -> Integer
getcigseqdiff a =
    let tcig = trimdcigar a
        cigmatchlen = sumSeqMatches $ mapcig tcig
        seqlen = B.length $ refseq a
    in  (fromIntegral seqlen) - cigmatchlen

getTrimdcigCoordDiff :: AlignedRead -> Integer
getTrimdcigCoordDiff a
    | tcig == "*" = 0
    | otherwise = diff
        where tcig = trimdcigar a
              tcmatchlen = sumRefMatches $ mapcig tcig
              trimdalen = (trimdendpos a) - (trimdpos a)
              diff =  trimdalen - tcmatchlen


-- 170206 filter nonmatching and zero-length CIGAR/sequence lengths
-- however: keep "*"
checkcigseqlen2 :: AlignedRead -> Bool
checkcigseqlen2 a
    | cigar a == "*" = True
    | (cigmatchlen == refseqlen) && (matchcnt > 0) = True
    | otherwise = False
        where matchcnt = sumMatches tcmap
              cigmatchlen = sum [ x | (x, y) <- tcmap, y == "M"
                                                    || y == "I"
                                                    || y == "S" ]
              tcmap = mapcig $ trimdcigar a
              refseqlen = genericLength $ B.unpack $ refseq a

-- 171018 new CIGAR sequence length check which only checks non-header AlignedRead recs
checkCigarSeqlen :: AlignedRead -> AlignedRead
checkCigarSeqlen a
    | (isheader a) = a
    | trimdcigarmap == Nothing = defaultAlignment
    | (cigseqlenHdrPassTest a) == True = a
    | otherwise = defaultAlignment
        where trimdcigarmap = safemapcig $ trimdcigar a

cigseqlenHdrPassTest :: AlignedRead -> Bool
cigseqlenHdrPassTest a
    -- | (isheader a) = True
    | (cigar a == "*") || (trimdcigar a == "*") = True -- 180222
    | (cigmatchlen == refseqlen) = True
    -- | (cigmatchlen == refseqlen) && (matchcnt > 0) = True
    | otherwise = False
        where matchcnt = sumMatches tcmap
              cigmatchlen = sum [ x | (x, y) <- tcmap, y == "M"
                                                    || y == "I"
                                                    || y == "S" ]
              tcmap = mapcig $ trimdcigar a
              refseqlen = genericLength $ B.unpack $ refseq a

-- translate CIGAR string into aligned interval with CIGAR op labels
-- 171018 safer version due to head on empty list error
safemapcig :: B.ByteString -> Maybe CigarMap
safemapcig cigstr
    | (length cigtups) > 0 = Just $ head cigtups
    | otherwise = Nothing
        where cigtups = U.rights $ [parseCigar cigstr] -- [(Integer, B.ByteString)]

-- translate CIGAR string into aligned interval with CIGAR op labels
mapcig :: B.ByteString -> CigarMap
mapcig cigstr =
    let cigtups = exrights $ parseCigar cigstr -- [(Integer, B.ByteString)]
    in cigtups

-- tuple map helper functions
mapfst :: Num a => (a -> c) -> (a, b) -> (c, b)
mapfst f (x, y) = (f x, y)

-- parsers and helper functions
readint = read :: String -> Integer
rdint = read :: String -> Int

parseint = A.parseOnly A.decimal

parsedbl = A.parseOnly A.double

parsesignedint i = floor $ exrights $ parsedbl i :: Integer

spaces = A.many1 A.skipSpace

--txtfieldp = A.takeTill A.isHorizontalSpace
txtfieldp = A.takeTill A.isSpace

optfieldsp = A.sepBy' txtfieldp A.space -- parses correctly

optfieldstotalp = A.takeTill (A.inClass "\r\n")

-- 171017 prevent entire SAM file being parsed into optional fields field of
-- first AlignedRead
optfieldpEOL = A.manyTill (A.satisfy (A.notInClass "\r\n")) A.endOfLine

optfieldsp2 = A.sepBy' optfieldpEOL A.space

optfieldsp3 = A.sepBy' txtfieldp A.space

optfieldp :: A.Parser B.ByteString
optfieldp = optnmp <|> optasp <|> optxsp

-- number of mismatches in alignment
optnmp :: A.Parser B.ByteString
optnmp = do
    id <- optidp
    skipcolon
    dtype <- itypep
    skipcolon
    cnt <- A.decimal
    let cntstr = B.pack $ show cnt
        nmoptstr = B.intercalate ":" [id, dtype, cntstr]
    return nmoptstr

-- alignment score of this alignment
optasp :: A.Parser B.ByteString
optasp = do
    id <- optidp
    skipcolon
    dtype <- itypep
    skipcolon
    score <- A.decimal
    let scorestr = B.pack $ show score
        asoptstr = B.intercalate ":" [id, dtype, scorestr]
    return asoptstr

-- alignment score of best alternate alignment (if AS is >= XS
-- then this aln is primary)
optxsp :: A.Parser B.ByteString
optxsp = do
    id <- optidp
    skipcolon
    dtype <- itypep
    skipcolon
    score <- A.decimal
    let scorestr = B.pack $ show score
        xsoptstr = B.intercalate ":" [id, dtype, scorestr]
    return xsoptstr

skipcolon = A.skipWhile (\x -> x == ':')

optidp :: A.Parser B.ByteString
optidp =    ("NM" >> return "NM")
        <|> ("AS" >> return "AS")
        <|> ("XS" >> return "XS")

itypep :: A.Parser B.ByteString
itypep = "i" >> return "i"

-- 171017 mID parsing (additional mID processing to be added to master as needed)
mIDp :: A.Parser B.ByteString
mIDp = do
    tag <- A.string "RX:Z"
    mid <- A.takeWhile (A.inClass "ACGTN")
    return mid

seqp :: A.Parser B.ByteString
seqp = A.takeWhile (A.inClass "ACTGatcgN")

parsemIDstring :: B.ByteString -> B.ByteString
parsemIDstring bs
    | length midvals <= 0 = "NOMIDSTRING"
    | otherwise = head midvals
        where midvals = U.rights $ (A.parseOnly mIDp) <$> (B.words bs)


intgr2int :: Integer -> Int
intgr2int n = fromIntegral n

-- convert decimal flag to binary representation for parsing flag
toBin :: Int -> String
toBin n = concatMap show $ digits 2 n

toBinInt :: String -> Int
toBinInt n = read n

-- possible data structure for SAM flag
data SAMFlag = SAMFlag { pairedRead :: Bool
                       , mateMapped :: Bool
                       , notMapped  :: Bool
                       , pairNotMapped :: Bool
                       , queryPlus :: Bool
                       , mateNeg :: Bool
                       , fstInPair :: Bool
                       , sndInPair :: Bool
                       , notPrimary :: Bool
                       , failedQC :: Bool
                       , dupRead :: Bool
                       , intflag :: Int
                       } deriving (Show, Eq)

readSAMFlag :: Int -> SAMFlag
readSAMFlag flag =
       SAMFlag { pairedRead = testBit flag 0
               , mateMapped = testBit flag 1
               , notMapped = testBit flag 2
               , pairNotMapped = testBit flag 3
               , queryPlus = testBit flag 4
               , mateNeg =  testBit flag 5
               , fstInPair = testBit flag 6
               , sndInPair = testBit flag 7
               , notPrimary = testBit flag 8
               , failedQC = testBit flag 9
               , dupRead = testBit flag 10
               , intflag = flag
               }

-- filter alignments by flag bit(s)
-- 180321
read1filter :: AlignedRead -> Bool
read1filter a = testBit (flag a) 6

primaryR1filter :: AlignedRead -> Bool
primaryR1filter a = (not $ testBit (flag a) 8)
                 && (read1filter a)

primaryR2filter :: AlignedRead -> Bool
primaryR2filter a =  (not $ testBit (flag a) 8)
                 && (not $ read1filter a)

-- 180322 filter PairedAln for records containing at least one alignment with
-- an intersection to at least one primer interval
collectPrimIntAlns :: [PairedAln] -> [PairedAln]
collectPrimIntAlns ps = filter anyPrimerIntAln ps

anyPrimerIntAln :: PairedAln -> Bool
anyPrimerIntAln p = (pintflag $ r1prim p)
                 || (pintflag $ r2prim p)
                 || (any pintflag $ r1secs p)
                 || (any pintflag $ r2secs p)


-- select element of nested vector
getcol :: Int -> V.Vector (V.Vector a) -> V.Vector a
getcol n txt = fmap (V.! n) txt

-- short notation for selecting a vector element by index
ix :: Int -> V.Vector a -> a
ix i v = (V.! i) v

-- parsing hack
exrights x = head $ U.rights $ [x]

-- safer exrights
exrights2 :: [a] -> Maybe a
exrights2 xs
    | length xs == 0 = Nothing
    | otherwise = Just (head xs)

getlengths seqs = fmap B.length seqs

------------------------------------------------------------------------------
----------------------  HIGH-LEVEL FUNCTIONS IN MAIN -------------------------
------------------------------------------------------------------------------
{--
--parseTrimSAM :: P.MonadResource m => B.ByteString -> P.ConduitM 
parseTrimSAM outfile = do
    pe <- CA.conduitParserEither hdralnparserEOL'
    let defltPosRang = CA.PositionRange (CA.Position 0 0 0) (CA.Position 0 0 0)
        (d, p) = U.fromRight (defltPosRang, defaultAlignment) pe
    -- B.writeFile outfile $ printAlignmentOrHdr p
    -- P.sinkList -- pass remainder of stream through function
    return p
--}

{--
-- 180329 try bundling parsing and trimming into one function
parseAndTrimPairSet :: CMap -> CMap -> B.ByteString -> [AlignedRead]
parseAndTrimPairSet fmp rmp bs =
    -- assumes header already parsed from the stream
    let pdaln = U.fromRight defaultPairedAln
                $ A.parseOnly parsePairedAlnsOrHdr bs
        trimd = trimprimerPairsE fmp rmp pdaln
    in sort $ (r1prim trimd)
            : (r2prim trimd)
            : ((r1secs trimd) ++ (r2secs trimd))
--}

-- 180409 expand PairedAln to [AlignedRead]
flattenPairedAln :: PairedAln -> [AlignedRead]
flattenPairedAln p = sort $ (r1prim p) : (r2prim p) : ((r1secs p) ++ (r2secs p))

-- 180409 add function to update RNEXT and PNEXT of alns w/ 0-trimmed
-- 180417 try removing any secondary alignments which are 0-trimmed
trimprimerPairsE :: CMap -> CMap -> PairedAln -> PairedAln
trimprimerPairsE fmap rmap p =
    let intpaln = addprimerintsPairedAln fmap rmap p
        updated = makeTrimmedUpdates $ trimPairedAlns intpaln
        nonprimZeroLengthRemoved = removeNonPrimaryZeroLengthAlignments updated
    in nonprimZeroLengthRemoved

-- element (single alignment)-wise primer trimming (for use with conduit)
trimprimersE :: CMap -> CMap -> AlignedRead -> AlignedRead
trimprimersE fmap rmap a =
    let intaln = addprimerints fmap rmap a
        trimd = trimalignment intaln
    in trimd

-- 171017 function to convert conduit ZipSink stream to output and write
-- output to file in SAM format.
printAlnStreamToFile :: P.MonadResource m => FilePath -> P.ConduitM AlignedRead c m ()
printAlnStreamToFile outfile = P.mapC printAlignmentOrHdr
                          P..| P.unlinesAsciiC
                          P..| P.sinkFile outfile

-- 180226 write RunStats to run log file (TODO: also print to stderr to allow
-- more flexible logging from caller of primerclip?)
-- NOTE: leaving out % mapped (until/unless we try to report change in %
-- mapped due to trimming alignments to zero length, resulting in change in
-- % mapped
writeRunStats :: FilePath -> RunStats -> IO ()
writeRunStats fp r = do
    let tot = B.pack $ show $ alnsTotal r
        mapd = B.pack $ show $ alnsMapped r
        trimd = B.pack $ show $ alnsTrimd r
        tt0 = B.pack $ show $ alnsTrimdToZero r
        tPct = B.pack $ show $ trimmedPct r
        mPct = B.pack $ show $ mappedPct r
        labels = [ "Total alignments processed:"
                 , "Total mapped alignments:"
                 , "Alignments trimmed by >= 1 base:"
                 , "Alignments trimmed to zero aligned length:"
                 , "% Alignments trimmed by >= 1 base:"
                 , "% Alignments mapped after trimming:"
                 ] :: [B.ByteString]
        statslines = zipWith (\l v -> B.concat [l,"\t",v])
                             labels
                             [tot, mapd, trimd, tt0, tPct, mPct]
        outbs = B.unlines statslines
        logfilename = genLogFilePath fp
    B.writeFile logfilename outbs

-- 180226 safely try to get name root of input SAM file to create log filename
genLogFilePath :: FilePath -> FilePath
genLogFilePath fp
    | (length parts) < 2 = "primerclip_runstats.log"
    | otherwise = nameroot
        where bs = B.pack fp
              parts = B.split '.' bs -- NOTE: should we have an alternate split char e.g. '.'?
              nameroot = B.unpack
                       $ B.append (head parts) "_primerclip_runstats.log"

-- 171017 calculate trim stats and print to stdout (TODO: print full stats to file)
calculateTrimStats :: P.ConduitM AlignedRead c (P.ResourceT IO) Integer
calculateTrimStats = P.filterC (\x -> trimdflag x) P..| P.lengthC

-- 180226 use ZipSink to calculate run stats and return a RunStats record for
-- printing to stdout, stderr, and to log file (purposely keeping the IO for
-- the RunStats outside Conduit pipeline for now)
calcRunStats = (calc <$> P.ZipSink P.lengthC
                     <*> P.ZipSink calcMappedCount
                     <*> P.ZipSink calculateTrimStats
                     <*> P.ZipSink calcTrimdToZero)
    where calc total mapd trimd trimd2z = RunStats total
                                                   mapd
                                                   trimd
                                                   trimd2z
                                                   ((fromIntegral trimd)
                                                    / (fromIntegral total)
                                                    * 100.0)
                                                   ((fromIntegral mapd)
                                                    / (fromIntegral total)
                                                    * 100.0)

calcMappedCount :: Integral i => P.ConduitM AlignedRead c (P.ResourceT IO) i
calcMappedCount = P.filterC (\x -> (mapped x) && (not $ trimdToZeroLength x))
                  P..| P.lengthC

calcTrimdToZero :: Integral i => P.ConduitM AlignedRead c (P.ResourceT IO) i
calcTrimdToZero = P.filterC (\x -> trimdToZeroLength x) P..| P.lengthC

-- 180226 calculate percentage of total alignments trimmed by >=1 bases
calcTrimmedPct :: P.Sink AlignedRead (P.ResourceT IO) Double
calcTrimmedPct = P.getZipSink (calc <$> P.ZipSink calculateTrimStats
                                    <*> P.ZipSink P.lengthC)
    where calc trimdcnt totalalns = (fromIntegral trimdcnt)
                                  / (fromIntegral totalalns) * 100.0

-- 180226 calculate percentage of alignments not mapped to reference after trimming
-- calcNotMappedPct :: P.Sink AlignedRead (P.ResourceT IO) Double
calcNotMappedPct = P.getZipSink (calc <$> P.ZipSink calcMappedCount
                                      <*> P.ZipSink P.lengthC)
    where calc mapcnt total = (fromIntegral mapcnt)
                            / (fromIntegral total) * 100.0


-- 180322 add primer intersections to AlignedRead as part of PairedAln record
-- NOTE: define instance of Functor for PairedAln type to define fmap??
addprimerintsPairedAln :: CMap -> CMap -> PairedAln -> PairedAln
addprimerintsPairedAln fpmap rpmap pa =
    let addpints = addprimerints fpmap rpmap
    in pa { r1prim = addpints $ r1prim pa
          , r2prim = addpints $ r2prim pa
          , r1secs = addpints <$> (r1secs pa)
          , r2secs = addpints <$> (r2secs pa)
          }

-- 170926 calculate and populate amplicon target BED field
addprimerints :: CMap -> CMap -> AlignedRead -> AlignedRead
addprimerints fpmap rpmap aln =
    let achr = rname aln
        start = pos aln
        end = endpos aln
        fbedmap = justchrmaps $ [M.lookup achr fpmap] -- possibly []
        rbedmap = justchrmaps $ [M.lookup achr rpmap] -- possibly []
        starthits = justbedmaps $ (bedmaplookup start) <$> fbedmap
        endhits = justbedmaps $ (bedmaplookup end) <$> rbedmap
        allhits = sort (starthits ++ endhits) -- [BedRecord]
        flagstatus = setpintflag allhits -- set True if alignment intersects one or more primers
    in aln { fint = starthits
           , rint = endhits
           , pintflag = flagstatus
           }

setpintflag hits
    | (length hits) > 0 = True
    | otherwise = False


---------- Functions to trim AlignedRead if it intersects primer(s) ----------
------------------------------------------------------------------------------

-- {--
-- 180322 update flag, rnext, and pnext field for each alignment in a PairedAln
-- record based on result of any primer trimming of that read's mate
updateTrimdPairFields :: PairedAln -> PairedAln
updateTrimdPairFields pa
    | bothPrimaryTrimd = updateR1nextfields $ updateR2nextfields pa
    | primaryR1trimd = updateR2nextfields pa
    | primaryR2trimd = updateR1nextfields pa
    | otherwise = pa
        where r1p = r1prim pa
              r2p = r2prim pa
              bothPrimaryTrimd = (trimdflag r1p)
                              && (trimdflag r2p)
                              -- && (mapped r1p)
                              -- && (mapped r2p)
              primaryR1trimd = (trimdflag r1p) -- && (mapped r1p)
              primaryR2trimd = (trimdflag r2p) -- && (mapped r2p)

-- called on PairedAln records where R2 primary alignment was trimmed
updateR1nextfields :: PairedAln -> PairedAln
updateR1nextfields pa =
    let trimdaln = r2prim pa
        trimdposR2 = trimdpos trimdaln
        trimdflagR2 = flag trimdaln
        nxtalns = (r1prim pa) : (r1secs pa)
        (newpr1:newsecr1s) = (\x -> x { pnext = trimdposR2 }) <$> nxtalns -- update pnext
    in pa { r1prim = newpr1, r1secs = newsecr1s }

updateR2nextfields :: PairedAln -> PairedAln
updateR2nextfields pa =
    let trimdaln = r1prim pa
        trimdposR1 = trimdpos trimdaln
        trimdflagR1 = flag trimdaln
        nxtalns = (r2prim pa) : (r2secs pa)
        (newpr2:newsecr2s) = (\x -> x { pnext = trimdposR1 }) <$> nxtalns -- update pnext
        -- NOTE: should MRNM be kept "*" for trimmed-to-0-length with mate still mapped???
    in pa { r2prim = newpr2, r2secs = newsecr2s }
--}

-- 180402 test whether paired alignment start is correctly updated for
-- primer-trimmed alignments
-- NOTE: use only on PairedAln w/ >=1 AlignedRead where trimdflag == True
validateTrimdPairAlnStart :: PairedAln -> Bool
validateTrimdPairAlnStart p = ((pnext $ r2prim p) == (trimdpos $ r1prim p))
                           && ((pnext $ r1prim p) == (trimdpos $ r2prim p))


-- 180322 trim each alignment in PairedAln record
trimPairedAlns :: PairedAln -> PairedAln
trimPairedAlns pa = PairedAln (trimalignment $ r1prim pa)
                              (trimalignment $ r2prim pa)
                              (trimalignment <$> (r1secs pa))
                              (trimalignment <$> (r2secs pa))

-- 180411 remove updateTrimdAlnFields and apply after MRNM resolution for 0-trimd alns
-- 180212 added addtrimtag function to append comment to optfields (CO:Z: SAM tag)
-- 180417 check alignment mapped before trying to trim (check mapped flag bit)
trimalignment :: AlignedRead -> AlignedRead
trimalignment a
    | (fint a /= []) && (rint a /= []) && mapdaln = btrimdaln
    | (fint a /= []) && (rint a == []) && mapdaln = ftrimdaln
    | (fint a == []) && (rint a /= []) && mapdaln = rtrimdaln
    | (fint a == []) && (rint a == []) = a { trimdcigar = cigar a }
    | otherwise = a { trimdcigar = cigar a }
         where ftrimdaln = trimfwd a
               rtrimdaln = trimrev a
               btrimdaln = trimboth a
               mapdaln = not $ flipTstBit 2 (flag a)

-- for alignments intersecting a "forward" primer only ( ref (+)-strand orientation)
trimfwd :: AlignedRead -> AlignedRead
trimfwd a =
    let newpos = bedend $ head $ fint a
        as = pos a
        fdiff = newpos - as -- handle alignment ending "inside" primer
        oldcigar = cigar a
        newcig = updateCigF fdiff oldcigar
        newcigmap = exrights $ parseCigar newcig
        tpos = (if (fdiff < 0) then as else newpos)
        zerolengthflag = boolZeroLengthCigar newcig -- 180411
        trimdaln = a { trimdpos = tpos
                     , trimdcigar = newcig
                     , trimdcigmap = newcigmap
                     , trimdflag = if (as /= tpos) then True else False
                     , trimdToZeroLength = zerolengthflag
                     }
    in trimdaln

-- for alignments intersecting a "reverse" primer only ( ref (+)-strand orientation)
trimrev :: AlignedRead -> AlignedRead
trimrev a =
    let newendpos = bedstart $ head $ rint a
        ae = endpos a
        rdiff = (ae - newendpos) -- handle alignment ending "inside" primer
        oldcigar = cigar a
        newcig = updateCigR rdiff oldcigar
        newcigmap = exrights $ parseCigar newcig
        tendpos = (if (rdiff < 0) then ae else newendpos)
        zerolengthflag = boolZeroLengthCigar newcig -- 180411
        trimdaln = a { trimdendpos = tendpos
                     , trimdcigar = newcig
                     , trimdcigmap = newcigmap
                     , trimdflag = if (ae /= tendpos) then True else False
                     , trimdToZeroLength = zerolengthflag
                     }
    in trimdaln

-- for alignments with primer intersections at both ends
trimboth :: AlignedRead -> AlignedRead
trimboth a =
    let newpos = bedend $ head $ fint a
        newendpos = bedstart $ head $ rint a
        as = pos a
        ae = endpos a
        fdiff = (newpos - as)
        rdiff = (ae - newendpos)
        oldcigar = cigar a
        newcig = updateCigB fdiff rdiff oldcigar
        newcigmap = exrights $ parseCigar newcig
        tpos = (if (fdiff < 0) then as else newpos)
        tendpos = (if (rdiff < 0) then ae else newendpos)
        zerolengthflag = boolZeroLengthCigar newcig -- 180411
        trimdaln = a { trimdpos = tpos
                     , trimdendpos = tendpos
                     , trimdcigar = newcig
                     , trimdcigmap = newcigmap
                     , trimdflag = if ((as /= tpos) || (ae /= tendpos))
                                   then True
                                   else False
                     , trimdToZeroLength = zerolengthflag
                     }
    in trimdaln

-- UPDATE 18-03-01 revert CIGAR trimming changes while debugging problems
updateCigF :: Integer -> B.ByteString -> B.ByteString
updateCigF fdiff cigar
    | snd (head cmap) == "*" = "*"
    | fdiffi <= 0 = cigar
    | ((nopadlen - fdiffi) > 0) = newcig -- 180219 DEBUG testing
    | otherwise = "*"
        where cmap = mapcig cigar
              grps = B.group $ expandcigar cmap
              nohardgrps = B.group $ expandcigar $ filter nohardclip cmap
              fHs = B.filter (== 'H') $ head grps -- maybe ""
              rHs = B.filter (== 'H') $ last grps -- maybe ""
              fSs = B.filter (== 'S') $ head nohardgrps
              rSs = B.filter (== 'S') $ last nohardgrps
              fdiffi = intgr2int fdiff
              cignoclip = filter noclip cmap
              cigexp = expandcigar2 cignoclip -- [(Int, B.ByteString)] assoc. list
              nopadlen = length cigexp
              (ftrimCigOps, remCigOps) = splitAt fdiffi cigexp -- 180223
              ftrimSlength = sumSoftClipCigOps ftrimCigOps
              newss = B.replicate ftrimSlength 'S'
              newcigarcore = B.append newss $ B.concat $ snd <$> remCigOps -- remcigraw
              newcigar = B.append (B.append fSs newcigarcore) rSs
              newfullcigar = B.append fHs (B.append newcigar rHs)
              newcig = contractcigar newfullcigar

updateCigR :: Integer -> B.ByteString -> B.ByteString
updateCigR rdiff cigar
    | snd (head cmap) == "*" = "*"
    | rdiffi <= 0 = cigar
    | ((nopadlen - rdiffi) > 0) = newcig -- 180219 DEBUGGING
    | otherwise = "*"
        where cmap = mapcig cigar
              grps = B.group $ expandcigar cmap
              nohardgrps = B.group $ expandcigar $ filter nohardclip cmap
              fHs = B.filter (== 'H') $ head grps -- maybe ""
              rHs = B.filter (== 'H') $ last grps -- maybe ""
              fSs = B.filter (== 'S') $ head nohardgrps
              rSs = B.filter (== 'S') $ last nohardgrps
              rdiffi = intgr2int rdiff
              cignoclip = filter noclip cmap
              cigexpR = expandcigar2 $ reverse cignoclip
              nopadlen = length cigexpR
              (rtrimCigOps, remCigOps) = splitAt rdiffi cigexpR -- 180223
              rtrimSlength = sumSoftClipCigOps rtrimCigOps
              newss = B.replicate rtrimSlength 'S'
              newcigarcore = B.append
                            (B.concat $ snd <$> (reverse remCigOps))
                             newss
              newcigar = B.append fSs (B.append newcigarcore rSs)
              newfullcigar = B.append fHs (B.append newcigar rHs)
              newcig = contractcigar newfullcigar

updateCigB :: Integer -> Integer -> B.ByteString -> B.ByteString
updateCigB fdiff rdiff cigar
    | snd (head cmap) == "*" = "*"
    | fdiffi <= 0 = updateCigR rdiff cigar
    | rdiffi <= 0 = updateCigF fdiff cigar
    | ((nopadlen - fdiffi - rdiffi) > 0) = newcig -- 180212 zero-match CIGAR strings fail picard validation
    | otherwise = "*"
        where cmap = mapcig cigar
              grps = B.group $ expandcigar cmap
              nohardgrps = B.group $ expandcigar $ filter nohardclip cmap
              fHs = B.filter (== 'H') $ head grps -- maybe ""
              rHs = B.filter (== 'H') $ last grps -- maybe ""
              fSs = B.filter (== 'S') $ head nohardgrps
              rSs = B.filter (== 'S') $ last nohardgrps
              fdiffi = intgr2int fdiff
              rdiffi = intgr2int rdiff
              -- 5p trim
              cignoclipf = filter noclip cmap
              cigexpf = expandcigar2 cignoclipf
              nopadlen = length cigexpf
              (ftrimCigOps, fremCigOps) = splitAt fdiffi cigexpf -- 180223
              ftrimSlength = sumSoftClipCigOps ftrimCigOps
              newfss = B.replicate ftrimSlength 'S'
              cigexpftrimd = reverse fremCigOps
              (rtrimCigOps, rremCigOps) = splitAt rdiffi cigexpftrimd -- 180223
              rtrimSlength = sumSoftClipCigOps rtrimCigOps
              newrss = B.replicate rtrimSlength 'S'
              totfss = B.append fSs newfss
              totrss = B.append rSs newrss
              newcigarcore = B.concat $ snd <$> (reverse rremCigOps)
              newcigar = B.append totfss (B.append newcigarcore totrss)
              newfullcigar = B.append fHs (B.append newcigar rHs)
              newcig = contractcigar newfullcigar

showcigar :: (Integer, B.ByteString) -> B.ByteString
showcigar cm =
    let cnt = B.pack $ show $ fst cm
    in B.append cnt (snd cm)

expandcigar :: CigarMap -> B.ByteString
expandcigar cmap = B.concat
                $ [ B.replicate (intgr2int n) $ B.head c
                  | (n, c) <- cmap ]

expandcigar2 :: CigarMap -> [(Integer, B.ByteString)]
expandcigar2 cmap =
    let ops = B.singleton <$> (B.unpack $ expandcigar cmap)
        opslen = genericLength ops
        ixs = [1..opslen] :: [Integer]
    in zipWith (,) ixs ops

-- convert expanded cigar string to standard CIGAR format (after trim update)
contractcigar :: B.ByteString -> B.ByteString
contractcigar longcig
    | longcig == "*" = "*"
    | otherwise = cigarstring
        where grpd = B.group longcig -- maintains order of groups
              shorten x = B.append (B.pack $ show $ B.length x)
                                   (B.singleton $ B.head x)
              cigarstring = B.concat $ shorten <$> grpd

-- 161023 sum all "S", "I", and "M" chars in trimmed sequence (drop all "D")
-- to get accurate CIGAR string
filtpadassoc :: [(Int, B.ByteString)] -> [(Int, B.ByteString)]
filtpadassoc asclist = filter nopadassoc asclist
    where nopadassoc = (\x -> (snd x) /= "D")

-- 161107 remove "D"s before making associative array ([(Int, B.ByteString)])
filtpadcmap :: CigarMap -> CigarMap
filtpadcmap cmap = filter nopad cmap
    where nopad = (\x -> (snd x) /= "D")

-- 161023 remove "H" and "D" from CigarMap (for trimming by soft clipping)
nopadding :: (Integer, B.ByteString) -> Bool
nopadding cp
    | (op /= "H") && (op /= "D") = True
    | otherwise = False
        where op = snd cp

-- 170201
removeDs :: [(Integer, B.ByteString)] -> [(Integer, B.ByteString)]
removeDs cs = filter (\x -> (snd x) /= "D") cs

countDs :: [(Integer, B.ByteString)] -> Integer
countDs cs = genericLength $ filter (\x -> (snd x) == "D") cs

nohardclip :: (Integer, B.ByteString) -> Bool
nohardclip op
    | (snd op) /= "H" = True
    | otherwise = False

noclip :: (Integer, B.ByteString) -> Bool
noclip op
    | ((snd op) /= "H") && ((snd op) /= "S") = True
    | otherwise = False

-- collect cigar string group (usually 'H' or 'S') into ordered list
getNclips :: [B.ByteString] -> Char -> [B.ByteString]
getNclips ciggrps cigchar = filter (\x -> B.length x > 0)
                          $ B.filter (== cigchar) <$> ciggrps

-- 171018 conduit-compatible filters for cigar string match count > 0
-- this function is designed to filter out alignments which have zero matches to
-- the reference after primer trimming
checknonzeroCigMatch :: AlignedRead -> Bool
checknonzeroCigMatch a
    | (isheader a) = True
    | (not $ mapped a) = True
    | (B.any (== 'M') (trimdcigar a)) = True
    | otherwise = False

-- 180213 check for "no real operator" (picard) in trimmed CIGAR string
-- 180222 several problems with invalid flag settings once CIGAR is cleared;
-- thoroughly test edge cases to ensure picard does not fail due to malformed
-- alignments
clearNonRealCigar :: AlignedRead -> AlignedRead
clearNonRealCigar a
    | (any (\x -> elem x ("MIDN" :: String)) (B.unpack $ trimdcigar a)) = a
    | otherwise = clearedCigAln
        where clearedCigAln = a { trimdcigar = "*"
                                -- , flag = zeroLenFlag
                                , trimdpos = 0
                                , trimdendpos = 0
                                -- , rnext = "*"
                                -- , pnext = -1
                                , tlen = 0
                                , mapqual = 0
                                , trimdToZeroLength = True
                                , mapped = False
                                }
              -- zeroLenFlag = setZeroLengthAlnFlag $ flag a


-- {--
-- 180409 clear flags only on trimmed-to-zero-length Alignments
-- (NOT their mapped pairs)
-- TODO: set mate MRNM of trimmed-to-zero-length alignments to ... ?
-- 180726 check case of one read not mapped pre-primerclip, and other read
-- trimmed to zero-length alignment
updateZeroTrimdPairFlags :: PairedAln -> PairedAln
updateZeroTrimdPairFlags pa
    | (r1zerotrimd && r2zerotrimd) = clearedBothNextmapflags
    | r1zerotrimd = clearedR2nextmapflags
    | r2zerotrimd = clearedR1nextmapflags
    | otherwise = pa
        where r1zerotrimd = (trimdToZeroLength r1p)
              r2zerotrimd = (trimdToZeroLength r2p)
              clearedBothNextmapflags = pa { r1prim = r1pZ
                                           , r2prim = r2pZ
                                           , r1secs = r1Zs
                                           , r2secs = r2Zs
                                           }
              clearedR1nextmapflags = pa { r1prim = r1pZ
                                         , r2prim = r2pMRNM
                                         , r1secs = r1Zs
                                         , r2secs = r2sMRNMs
                                         }
              clearedR2nextmapflags = pa { r1prim = r1pMRNM
                                         , r2prim = r2pZ
                                         , r1secs = r1sMRNMs
                                         , r2secs = r2Zs
                                         }
              (r1p, r2p, r1s, r2s) = ( (r1prim pa)
                                     , (r2prim pa)
                                     , (r1secs pa)
                                     , (r2secs pa) )
              clrFlagMapBits x = x { flag = setZeroLengthPairFlag (flag x) }
              r1pZ = clrFlagMapBits r1pMRNM
              r2pZ = clrFlagMapBits r2pMRNM
              r1Zs = clrFlagMapBits <$> r1sMRNMs
              r2Zs = clrFlagMapBits <$> r2sMRNMs
              r1pMRNM = r1prim newMRNMp
              r2pMRNM = r2prim newMRNMp
              r1sMRNMs = r1secs newMRNMp
              r2sMRNMs = r2secs newMRNMp
              newMRNMp = setMateRname pa
--}

-- 180416 setMateRname must only update RNAME when mate was mapped in input SAM
setMateRname :: PairedAln -> PairedAln
setMateRname p
    | (r1mateorigmapped && r2mateorigmapped) = newMRNMprdaln
    | r1mateorigmapped = newr1MRNM
    | r2mateorigmapped = newr2MRNM
    | otherwise = p
        where r1mateorigmapped = not $ flipTstBit 2 (flag r2p)
              r2mateorigmapped = not $ flipTstBit 2 (flag r1p)
              newr1MRNM = p { r1prim = newr1p, r1secs = newr1secs }
              newr2MRNM = p { r2prim = newr2p, r2secs = newr2secs }
              newMRNMprdaln = p { r1prim = newr1p
                                , r2prim = newr2p
                                , r1secs = newr1secs
                                , r2secs = newr2secs
                                }
              newr1p = setRname r1p r2p
              newr2p = setRname r2p r1p
              newr1secs = (flip setRname r2p) <$> r1s
              newr2secs = (flip setRname r1p) <$> r2s
              setRname r m = r { rnext = (B.pack $ show $ rname m) }
              (r1p, r2p, r1s, r2s) = ( (r1prim p)
                                     , (r2prim p)
                                     , (r1secs p)
                                     , (r2secs p) )

-- {--
-- 180416 clear RNEXT and PNEXT for (primary) pair if primary trimmed to 0-length
updateZeroTrimdPairFields :: PairedAln -> PairedAln
updateZeroTrimdPairFields p
    | (primR1zerotrimmed && primR2zerotrimmed) = clearR2primNextFields
                                               $ clearR1primNextFields p
    | primR1zerotrimmed = clearR2primNextFields p
    | primR2zerotrimmed = clearR1primNextFields p
    | otherwise = p
        where primR1zerotrimmed = trimdToZeroLength $ r1prim p
              primR2zerotrimmed = trimdToZeroLength $ r2prim p
--}

-- {--
-- 180726 IFF (mate read not mapped by bwa && read is trimmed to 0-length)
-- then clear RNAME, POS, RNEXT (MRNM), and PNEXT for both reads
clearR1primNextFields :: PairedAln -> PairedAln
clearR1primNextFields p
    | not $ mapped $ r1prim p = clearNamesAndPositions p
    | otherwise = clearedNextFields
        where r1p = r1prim p
              newr1prim = r1p { rnext = "*"
                              , pnext = 0
                              , tlen = 0
                              }
              clearedNextFields = p { r1prim = newr1prim }

clearR2primNextFields :: PairedAln -> PairedAln
clearR2primNextFields p
    | not $ mapped $ r2prim p = clearNamesAndPositions p
    | otherwise = clearedNextFields
        where r2p = r2prim p
              newr2prim = r2p { rnext = "*"
                              , pnext = 0
                              , tlen = 0
                              }
              clearedNextFields = p { r2prim = newr2prim }

clearNamesAndPositions :: PairedAln -> PairedAln
clearNamesAndPositions p =
    let r1p = r1prim p
        r2p = r2prim p
        newr1p = r1p { rname = NONE
                     , pos = 0
                     , rnext = "*"
                     , pnext = 0
                     , tlen = 0
                     }
        newr2p = r2p { rname = NONE
                     , pos = 0
                     , rnext = "*"
                     , pnext = 0
                     , tlen = 0
                     }
    in p { r1prim = newr1p, r2prim = newr2p }

--}

{--
clearR1primNextFields :: PairedAln -> PairedAln
clearR1primNextFields p =
    let r1p = r1prim p
        newr1prim = r1p { rnext = "*"
                        , pnext = 0
                        , tlen = 0
                        }
    in p { r1prim = newr1prim }

clearR2primNextFields :: PairedAln -> PairedAln
clearR2primNextFields p =
    let r2p = r2prim p
        newr2prim = r2p { rnext = "*"
                        , pnext = 0
                        , tlen = 0
                        }
    in p { r2prim = newr2prim }
--}

-- 180423 add setProperInsertSizeRange step (TESTING)
-- modifications
-- TODO: consolidate ad-hoc updates to trimmed alignments,
--       and make setProperInsertSizeRange and range limits cmd line options
makeTrimmedUpdates :: PairedAln -> PairedAln
makeTrimmedUpdates pa = setProperInsertSizeRange (-1200) (1200)
                      $ updatePairedAlnTrimdFields
                      $ updateZeroTrimdPairFields
                      $ updateZeroTrimdPairFlags
                      $ updateTrimdPairFields pa

-- {--
makeMRNMexplicit :: PairedAln -> PairedAln
makeMRNMexplicit p
    | r1zerotrimdR2mapped || r2zerotrimdR1mapped = explicitMRNM
    | otherwise = p
        where r1zerotrimdR2mapped = (trimdToZeroLength $ r1prim p)
                                 || (any (\x -> trimdToZeroLength x) (r1secs p))
              r2zerotrimdR1mapped = (trimdToZeroLength $ r2prim p)
                                 || (any (\x -> trimdToZeroLength x) (r2secs p))
              explicitMRNM = p { r1prim = newr1p
                               , r1secs = newr1s
                               , r2prim = newr2p
                               , r2secs = newr2s
                               }
              -- explicitR2MRNM = p { r2prim = newr2p, r2secs = newr2s }
              newr1p = r1p { rnext = r2pRNAME }
              newr2p = r2p { rnext = r1pRNAME }
              newr1s = (\x -> x { rnext = r2pRNAME }) <$> r1s
              newr2s = (\x -> x { rnext = r1pRNAME }) <$> r2s
              r1p = r1prim p
              r2p = r2prim p
              r1s = r1secs p
              r2s = r2secs p
              r1pRNAME = B.pack $ show $ rname $ r1prim p
              r2pRNAME = B.pack $ show $ rname $ r2prim p
--}

-- 180320 clear supp. alignment bit
setZeroLengthAlnFlag :: Int -> Int
setZeroLengthAlnFlag flg
    | flipTstBit 0 flg = pairedZeroLengthFlag
    | otherwise = nopairZeroLengthFlag
        where pairedZeroLengthFlag = flipClrBit 11
                                   $ flipClrBit 8
                                   $ flipSetBit 2
                                   $ flipClrBit 1 flg
              nopairZeroLengthFlag = flipClrBit 11
                                   $ flipClrBit 8
                                   $ flipSetBit 2
                                   $ flipClrBit 1 flg

-- 180409 set mate not mapped bit
-- 180417 clear bit 1 (read mapped in proper pair)
setZeroLengthPairFlag :: Int -> Int
setZeroLengthPairFlag flg = flipSetBit 3
                           $ flipClrBit 1 flg

-- 180423 set/confirm bit 1 (read mapped in proper pair)
setProperPairMapFlagBit :: Int -> Int
setProperPairMapFlagBit flg = flipSetBit 1 flg

-- 180411 test for zero-trimmed alignment based on trimmed CIGAR
-- NOTE: prior to converting zero-match CIGAR to "*"
boolZeroLengthCigar :: B.ByteString -> Bool
boolZeroLengthCigar cigar = not
    $ any (\x -> elem x ("MIDN" :: String)) (B.unpack cigar)

-- 180411 apply updateTrimdAlnFields to all alignments in PairedAln
updatePairedAlnTrimdFields :: PairedAln -> PairedAln
updatePairedAlnTrimdFields p = PairedAln (updateTrimdAlnFields $ r1prim p)
                                         (updateTrimdAlnFields $ r2prim p)
                                         (updateTrimdAlnFields <$> r1secs p)
                                         (updateTrimdAlnFields <$> r2secs p)

updateTrimdAlnFields :: AlignedRead -> AlignedRead
updateTrimdAlnFields a
    | (any (\x -> elem x ("MIDN" :: String)) (B.unpack $ trimdcigar a))
        = trimdAln
    | (trimdflag a) = trimdToZero
    | otherwise = a -- no clipping due to no primer intersections for alignment
        where trimdAlnTag = "CO:Z:primer_trimmed" :: B.ByteString
              trimdToZeroTag = "CO:Z:zero_alignment_length_after_primer_trim" :: B.ByteString
              trimdAln = a { optfields = B.concat [ (optfields a)
                                                  , "\t", trimdAlnTag
                                                  ] }
              trimdToZero = a { optfields = B.concat [ (optfields a)
                                                     , "\t", trimdToZeroTag
                                                     ]
                              , rname = NONE
                              , trimdcigar = "*"
                              , trimdpos = -1
                              , trimdendpos = 0
                              , mapqual = 0
                              , mapped = False
                              , flag = setZeroLengthAlnFlag (flag a) -- 180409
                              }

-- 180417 remove any secondary alignments which are trimmed to zero-length
-- TODO: handle this inside makeTrimmedUpdates once we know we're keeping this step
removeNonPrimaryZeroLengthAlignments :: PairedAln -> PairedAln
removeNonPrimaryZeroLengthAlignments p =
    let (r1s, r2s) = ((r1secs p), (r2secs p))
        newr1s = filter (\x -> not $ trimdToZeroLength x) r1s
        newr2s = filter (\x -> not $ trimdToZeroLength x) r2s
    in p { r1secs = newr1s, r2secs = newr2s }

-- 180423 ensure flag bit 1 is set if trimmed PairedAln aligned insert size
-- is within specific range
-- TODO: make command line switch to enable/disable this option, and an arg
-- to specify the insert size range to use for keeping bit 1 set
-- NOTE: this function checks that bit 1 (read mapped in proper pair) is set
-- for PairedAln where both reads are mapped and insert size range is between
-- min and max (to ensure super-amplicon reads are included in variant calling).
setProperInsertSizeRange :: Integer -> Integer -> PairedAln -> PairedAln
setProperInsertSizeRange minsz maxsz p
    | inrange && pairmapped = ensureBit1Set
    | otherwise = p
        where ensureBit1Set = setMapdProperPairBit p
              inrange = checkInsertSize minsz maxsz p
              pairmapped = ((mapped $ r1prim p) && (mapped $ r2prim p))
                        -- || ((mapped $ r1prim p) && (any mapped (r2secs p))
                        -- || ((any mapped (r1secs p)) && (mapped $ r2prim p)

-- we're assuming that both R1 and R2 primary alignments are mapped after trimming.
-- TODO: test whether, after trimming, there are pairs where a secondary alignment
-- is mapped, but the primary alignment for that read is not mapped (presumably
-- due to trimming-to-zero-length...)
setMapdProperPairBit :: PairedAln -> PairedAln
setMapdProperPairBit p =
    let (r1p, r2p, _, _) = pairedAlnToTuple p
        newr1p = r1p { flag = (setProperPairMapFlagBit $ flag r1p) }
        newr2p = r2p { flag = (setProperPairMapFlagBit $ flag r2p) }
    in p { r1prim = newr1p, r2prim = newr2p }

checkInsertSize :: Integer -> Integer -> PairedAln -> Bool
checkInsertSize minsz maxsz p
    | (pairedmin >= minsz) && (pairedmax <= maxsz) = True
    | otherwise = False
        where pairedmin = minimum $ tlen <$> alns
              pairedmax = maximum $ tlen <$> alns
              alns = (r1prim p) : (r2prim p) : (join [(r1secs p), (r2secs p)])

-- flip setBit and clearBit args for clearer syntax
flipSetBit = flip setBit
flipClrBit = flip clearBit
flipTstBit = flip testBit

-- 180212 append CO:Z tag indicating alignment was trimmed by >= 1 base, and
-- also a warning if trimming removed all non-clipped bases from alignment
-- (alignment == primer sequence)
-- NOTE: this function also sets bit 2 and clears bit 1 of the SAM FLAG if
-- trimming results in zero-length alignment
addtrimtag :: AlignedRead -> AlignedRead
addtrimtag a
    | (trimdflag a) && ((trimdcigar a) /= "*") = trimdAln
    | (trimdflag a) = trimdToZero
    | otherwise = a -- no clipping due to no primer intersections for alignment
        where trimdAlnTag = "CO:Z:primer_trimmed" :: B.ByteString
              trimdToZeroTag = "CO:Z:zero_alignment_length_after_primer_trim" :: B.ByteString
              trimdAln = a { optfields = B.concat [ (optfields a)
                                                  , "\t", trimdAlnTag
                                                  ] }
              trimdToZero = a { optfields = B.concat [ (optfields a)
                                                     , "\t", trimdToZeroTag
                                                     ]
                              , flag = zeroLenFlag
                              }
              zeroLenFlag = setBit (clearBit (flag a) 1) 2 -- no longer mapped

-- 161017 create I.IntMap Int B.ByteString from (Int, BedRecord)
-- NOTE: this is an attempt at fast lookup of primer intervals
--       key represent chromosome coordinate (IntMap for each chromosome)
makechrbedmap :: BED -> M.Map UChr (I.IntMap BedRecord)
makechrbedmap bs =
    let blist = V.toList bs
        bbychr = groupBy (\x y -> bedchr x == bedchr y) blist -- [[BedRecord]]
        bchrs = bedchr <$> (head <$> bbychr) -- [UChr]
        chrintmaps = makebedmap <$> bbychr -- [I.IntMap BedRecord]
        tups = zipWith (,) bchrs chrintmaps
    in M.fromList tups

makebedmap :: [BedRecord] -> I.IntMap BedRecord
makebedmap bs =
    let tups = join $ makeprimertups <$> bs
    in I.fromList tups

-- expand primer BED interval
-- TODO: make target BED interval padding configurable by command line option
makeprimertups :: BedRecord -> [(Int, BedRecord)]
makeprimertups b =
    let ps = (bedstart b) - 2 -- pad border with arbitrary 2 bases for now
        pe = (bedend b) + 2 -- pad border for now
        rng = intgr2int $ (pe - ps) + 1
        recs = replicate rng $ b
        ints = intgr2int <$> [ps..pe]
    in zipWith (,) ints recs

-- IntMap lookup with default result when key not found (represents no primer
-- at that position)
bedmaplookup :: Integer -> I.IntMap BedRecord -> Maybe BedRecord
bedmaplookup pos bmap = I.lookup (intgr2int pos) bmap

justbedmaps :: [Maybe BedRecord] -> [BedRecord]
justbedmaps ms = catMaybes ms

justchrmaps :: [Maybe (I.IntMap BedRecord)] -> [I.IntMap BedRecord]
justchrmaps mcmaps = catMaybes mcmaps

createprimerbedmaps :: Opts -> IO ( M.Map UChr (I.IntMap BedRecord)
                                  , M.Map UChr (I.IntMap BedRecord) )
createprimerbedmaps args = case (bedpeformat args) of
    False -> do
        m <- getMasterFile $ incoordsfile args
        let fm = makechrbedmap $ masterToFPrimerBED m
            rm = makechrbedmap $ masterToRPrimerBED m
        return (fm, rm)
    True  -> do
        bedpe <- readBEDPE $ incoordsfile args
        let fb = bedpeToFbed <$> bedpe
            rb = bedpeToRbed <$> bedpe
            fmpFromBedPE = makechrbedmap $ V.fromList fb
            rmpFromBedPE = makechrbedmap $ V.fromList rb
        return (fmpFromBedPE, rmpFromBedPE)
    -- _        -> error "(!) Incorrect primer coord input string"

-- 180206 split BEDPE to forward and reverse BedRecords
bedpeToFbed :: BEDPE -> BedRecord
bedpeToFbed b = BedRecord (chr1 b) (start1 b) (end1 b)
                          (B.append (bedpename b) "_F")

bedpeToRbed :: BEDPE -> BedRecord
bedpeToRbed b = BedRecord (chr2 b) (start2 b) (end2 b)
                          (B.append (bedpename b) "_R")

-- 171017 convert AlignedRead to output string format, handling either header
-- information stored in an AlignedRead record or an actual alignment
-- NOTE: this simplifies parsing the input SAM file using Conduit
printAlignmentOrHdr :: AlignedRead -> B.ByteString
printAlignmentOrHdr a
    | (isheader a) = headerstr
    | otherwise = alnstr
        where rn = B.pack $ show $ rname a
              tpos = if (rn == "*") then 0 else (trimdpos a) + 1 -- SAM format
              q = qname a
              f = B.pack $ show $ flag a
              p = B.pack $ show $ tpos
              mq = B.pack $ show $ mapqual a
              c = trimdcigar a -- 170202 prints trimmed alignment CIGAR
              rnxt = rnext a
              pnxt = if (rnxt == "*") then 0 else (pnext a) + 1 -- index shift handling at 0
              pnxtBS = B.pack $ show $ pnxt -- 180213
              tl = B.pack $ show $ tlen a
              sq = refseq a
              bq = basequal a
              optfs = optfields a
              alnstr = B.intercalate "\t" [ q, f, rn, p, mq, c, rnxt,
                                            pnxtBS, tl, sq, bq, optfs ]
              headerstr = B.reverse $ B.drop 1 $ B.reverse
                        $ B.unlines $ headerstrings a

printAlignment :: AlignedRead -> B.ByteString
printAlignment a =
    let rn = B.pack $ show $ rname a
        tpos = if (rn == "*") then 0 else (trimdpos a) + 1 -- SAM format rules
        q = qname a
        f = B.pack $ show $ flag a
        p = B.pack $ show $ tpos
        mq = B.pack $ show $ mapqual a
        c = trimdcigar a -- 170202 prints trimmed alignment CIGAR
        rnxt = rnext a
        pnxt = if (rnxt == "*") then 0 else (pnext a) + 1 -- index shift handling at 0
        pnxtBS = B.pack $ show $ pnxt -- 180213
        tl = B.pack $ show $ tlen a
        sq = refseq a
        bq = basequal a
        optfs = optfields a
        alnstr = B.intercalate "\t" [ q, f, rn, p, mq, c, rnxt,
                                      pnxtBS, tl, sq, bq, optfs ]
    in alnstr

-- max(0,i)
checkpos :: Integer -> Integer
checkpos i
    | i < 0 = 0
    | otherwise = i

-- 170201 new trimming arithmetic adjusts read coords to reference
adjustcrds :: [(Integer, B.ByteString)] -> [(Integer, B.ByteString)]
adjustcrds cigs = scanl1 shiftcrds cigs

-- 'H' and 'S' should be removed from CIGAR string before calling this function
-- this function adjusts read positions to align with reference (relative)
-- coordinates, allowing correct calculation of the primer-trimmed CIGAR string
shiftcrds :: (Integer, B.ByteString) -> (Integer, B.ByteString)
        -> (Integer, B.ByteString)
shiftcrds a b
    | thiscig == "I" = adjustposI
    | thiscig == "D" = adjustposD
    | otherwise = keeppos
    where lastpos = fst a
          thispos = fst b
          thiscig = snd b
          adjustposI = (lastpos, thiscig)
          adjustposD = (lastpos + 2, thiscig)
          keeppos = (lastpos + 1, thiscig)

taketrim :: Integer -> [(Integer, B.ByteString)] -> [(Integer, B.ByteString)]
taketrim cnt cs = takeWhile (\x -> (fst x) <= cnt) cs

-- 170203 try to fix edge case where D is first base in remaining CIGAR string
trimrem :: Integer -> [(Integer, B.ByteString)] -> [(Integer, B.ByteString)]
trimrem cnt cs
    | length rem == 0 = []
    | (snd $ head rem) == "D" = tail rem
    | otherwise = rem
        where rem = dropWhile (\x -> (fst x) <= cnt) cs

-- 170510 print warning if chromosome name formats from SAM and master input
-- files do not match
checkChromNameMatchStatus :: Header -> BED -> IO ()
checkChromNameMatchStatus hdr bed = do
    let bedchrs = S.toList $ S.fromList $ bedchr <$> (V.toList bed) -- set of BED chrs
        hdrchrs = S.toList $ S.fromList $ U.rights $ A.parseOnly hdrchromp
                <$> ((hdr !!) <$> [1..25]) -- set of SAM hdr chrs
        matches = filter (\x -> x /= NONE) $ intersect bedchrs hdrchrs
    if (length matches) >= 1
        then putStrLn "SAM and BED chromosome name formats match."
        else error "ERROR: different chromosome name formats in SAM and BED!"

-- 180223 find aln by qname
findByQname :: B.ByteString -> [AlignedRead] -> [AlignedRead]
findByQname name as = filter (\x -> (qname x) == name) as

getNonHeaderAlns :: [AlignedRead] -> [AlignedRead]
getNonHeaderAlns as = filter (not . isheader) as

-- 180423 PairedAln to (r1primary, r2primary, r1seconds, r2seconds) 4-tuple
pairedAlnToTuple :: PairedAln ->
    (AlignedRead, AlignedRead, [AlignedRead], [AlignedRead])
pairedAlnToTuple p = ((r1prim p), (r2prim p), (r1secs p), (r2secs p))

-- end of Library

{--
-- UPDATE 18-03-01 revert CIGAR trimming changes while debugging problems
updateCigF :: Integer -> B.ByteString -> B.ByteString
updateCigF fdiff cigar
    | snd (head cmap) == "*" = "*"
    | fdiffi <= 0 = cigar
    | ((nopadlen - fdiffi) == 0) = "*" -- 180320
    | ((nopadlen - fdiffi) > 0) = newcig -- 180320
    | otherwise = "*"
        where cmap = mapcig cigar
              grps = B.group $ expandcigar cmap
              nohardgrps = B.group $ expandcigar $ filter nohardclip cmap
              fHs = B.filter (== 'H') $ head grps -- maybe ""
              rHs = B.filter (== 'H') $ last grps -- maybe ""
              fSs = B.filter (== 'S') $ head nohardgrps
              rSs = B.filter (== 'S') $ last nohardgrps
              fdiffi = intgr2int fdiff
              cignoclip = filter noclip cmap
              cigexp = expandcigar2 cignoclip -- [(Int, B.ByteString)] assoc. list
              nopadlen = length cigexp
              adjcig = adjustcrds cigexp
              ftrim = taketrim fdiff adjcig
              trimDs = countDs ftrim
              trimDcorr = intgr2int trimDs
              adjix = (genericLength ftrim) + trimDs
              remcigraw = drop trimDcorr $ trimrem adjix cigexp
              newss = B.replicate (length ftrim) 'S'
              newcigarcore = B.append newss $ B.concat $ snd <$> remcigraw
              newcigar = B.append (B.append fSs newcigarcore) rSs
              newfullcigar = B.append fHs (B.append newcigar rHs)
              newcig = contractcigar newfullcigar

updateCigR :: Integer -> B.ByteString -> B.ByteString
updateCigR rdiff cigar
    | snd (head cmap) == "*" = "*"
    | rdiffi <= 0 = cigar
    | ((nopadlen - rdiffi) == 0) = "*" -- 180320 DEBUGGING
    | ((nopadlen - rdiffi) > 0) = newcig -- 180320 DEBUGGING
    | otherwise = "*" -- 180320 allow all 'S' trimmed CIGAR (diff == nopadlen)
        where cmap = mapcig cigar
              grps = B.group $ expandcigar cmap
              nohardgrps = B.group $ expandcigar $ filter nohardclip cmap
              fHs = B.filter (== 'H') $ head grps -- maybe ""
              rHs = B.filter (== 'H') $ last grps -- maybe ""
              fSs = B.filter (== 'S') $ head nohardgrps
              rSs = B.filter (== 'S') $ last nohardgrps
              rdiffi = intgr2int rdiff
              cignoclip = filter noclip cmap
              cigexpR = expandcigar2 $ reverse cignoclip
              nopadlen = length cigexpR
              adjcigR = adjustcrds cigexpR
              rtrim = taketrim rdiff adjcigR
              trimDs = countDs rtrim
              trimDcorr = intgr2int trimDs
              adjix = (genericLength rtrim) + trimDs
              remcigraw = drop trimDcorr $ trimrem adjix cigexpR
              remDs = intgr2int $ countDs remcigraw
              remcigDadj = reverse $ drop remDs remcigraw
              newss = B.replicate ((length rtrim)
                                  + trimDcorr
                                  + remDs) 'S'
              newcigarcore = B.append (B.concat $ snd <$> remcigDadj) newss
              newcigar = B.append fSs (B.append newcigarcore rSs)
              newfullcigar = B.append fHs (B.append newcigar rHs)
              newcig = contractcigar newfullcigar

updateCigB :: Integer -> Integer -> B.ByteString -> B.ByteString
updateCigB fdiff rdiff cigar
    | snd (head cmap) == "*" = "*"
    | fdiffi <= 0 = updateCigR rdiff cigar
    | rdiffi <= 0 = updateCigF fdiff cigar
    | ((nopadlen - fdiffi - rdiffi) == 0) = "*" -- 180320
    | ((nopadlen - fdiffi - rdiffi) > 0) = newcig -- 180320
    | otherwise = "*"
        where cmap = mapcig cigar
              grps = B.group $ expandcigar cmap
              nohardgrps = B.group $ expandcigar $ filter nohardclip cmap
              fHs = B.filter (== 'H') $ head grps -- maybe ""
              rHs = B.filter (== 'H') $ last grps -- maybe ""
              fSs = B.filter (== 'S') $ head nohardgrps
              rSs = B.filter (== 'S') $ last nohardgrps
              fdiffi = intgr2int fdiff
              rdiffi = intgr2int rdiff
              -- 5p trim
              cignoclipf = filter noclip cmap
              cigexpf = expandcigar2 cignoclipf
              nopadlen = length cigexpf
              adjcigf = adjustcrds cigexpf
              ftrim = taketrim fdiff adjcigf
              ftrimDs = countDs ftrim
              ftrimDcorr = intgr2int ftrimDs
              adjixf = (genericLength ftrim) + ftrimDs
              remcigrawf = drop ftrimDcorr $ trimrem adjixf cigexpf
              newfss = B.replicate (length ftrim) 'S'
              -- 3p trim3p
              cigexpr = zipWith (\x (i, j) -> (x, j))
                                [1..]
                                (reverse remcigrawf)
              adjcigr = adjustcrds cigexpr
              rtrim = taketrim rdiff adjcigr
              rtrimDs = countDs rtrim
              rtrimDcorr = intgr2int rtrimDs
              adjixr = (genericLength rtrim) + rtrimDs
              remcigrawr = drop rtrimDcorr $ trimrem adjixr cigexpr
              remDsr = intgr2int $ countDs remcigrawr
              remcigrDadj = reverse $ drop remDsr remcigrawr
              newrss = B.replicate
                        ((length rtrim) + rtrimDcorr + remDsr + ftrimDcorr) 'S'
              totfss = B.append fSs newfss
              totrss = B.append rSs newrss
              newcigarcore = B.concat $ snd <$> remcigrDadj
              newcigar = B.append totfss (B.append newcigarcore totrss)
              newfullcigar = B.append fHs (B.append newcigar rHs)
              newcig = contractcigar newfullcigar
--}

{--
    chr1_KI270706v1_random
    chr1_KI270707v1_random
    chr1_KI270708v1_random
    chr1_KI270709v1_random
    chr1_KI270710v1_random
    chr1_KI270711v1_random
    chr1_KI270712v1_random
    chr1_KI270713v1_random
    chr1_KI270714v1_random
    chr2_KI270715v1_random
    chr2_KI270716v1_random
    chr3_GL000221v1_random
    chr4_GL000008v2_random
    chr5_GL000208v1_random
    chr9_KI270717v1_random
    chr9_KI270718v1_random
    chr9_KI270719v1_random
    chr9_KI270720v1_random
    chr11_KI270721v1_random
    chr14_GL000009v2_random
    chr14_GL000225v1_random
    chr14_KI270722v1_random
    chr14_GL000194v1_random
    chr14_KI270723v1_random
    chr14_KI270724v1_random
    chr14_KI270725v1_random
    chr14_KI270726v1_random
    chr15_KI270727v1_random
    chr16_KI270728v1_random
    chr17_GL000205v2_random
    chr17_KI270729v1_random
    chr17_KI270730v1_random
    chr22_KI270731v1_random
    chr22_KI270732v1_random
    chr22_KI270733v1_random
    chr22_KI270734v1_random
    chr22_KI270735v1_random
    chr22_KI270736v1_random
    chr22_KI270737v1_random
    chr22_KI270738v1_random
    chr22_KI270739v1_random
    chrY_KI270740v1_random
    chrUn_KI270302v1
    chrUn_KI270304v1
    chrUn_KI270303v1
    chrUn_KI270305v1
    chrUn_KI270322v1
    chrUn_KI270320v1
    chrUn_KI270310v1
    chrUn_KI270316v1
    chrUn_KI270315v1
    chrUn_KI270312v1
    chrUn_KI270311v1
    chrUn_KI270317v1
    chrUn_KI270412v1
    chrUn_KI270411v1
    chrUn_KI270414v1
    chrUn_KI270419v1
    chrUn_KI270418v1
    chrUn_KI270420v1
    chrUn_KI270424v1
    chrUn_KI270417v1
    chrUn_KI270422v1
    chrUn_KI270423v1
    chrUn_KI270425v1
    chrUn_KI270429v1
    chrUn_KI270442v1
    chrUn_KI270466v1
    chrUn_KI270465v1
    chrUn_KI270467v1
    chrUn_KI270435v1
    chrUn_KI270438v1
    chrUn_KI270468v1
    chrUn_KI270510v1
    chrUn_KI270509v1
    chrUn_KI270518v1
    chrUn_KI270508v1
    chrUn_KI270516v1
    chrUn_KI270512v1
    chrUn_KI270519v1
    chrUn_KI270522v1
    chrUn_KI270511v1
    chrUn_KI270515v1
    chrUn_KI270507v1
    chrUn_KI270517v1
    chrUn_KI270529v1
    chrUn_KI270528v1
    chrUn_KI270530v1
    chrUn_KI270539v1
    chrUn_KI270538v1
    chrUn_KI270544v1
    chrUn_KI270548v1
    chrUn_KI270583v1
    chrUn_KI270587v1
    chrUn_KI270580v1
    chrUn_KI270581v1
    chrUn_KI270579v1
    chrUn_KI270589v1
    chrUn_KI270590v1
    chrUn_KI270584v1
    chrUn_KI270582v1
    chrUn_KI270588v1
    chrUn_KI270593v1
    chrUn_KI270591v1
    chrUn_KI270330v1
    chrUn_KI270329v1
    chrUn_KI270334v1
    chrUn_KI270333v1
    chrUn_KI270335v1
    chrUn_KI270338v1
    chrUn_KI270340v1
    chrUn_KI270336v1
    chrUn_KI270337v1
    chrUn_KI270363v1
    chrUn_KI270364v1
    chrUn_KI270362v1
    chrUn_KI270366v1
    chrUn_KI270378v1
    chrUn_KI270379v1
    chrUn_KI270389v1
    chrUn_KI270390v1
    chrUn_KI270387v1
    chrUn_KI270395v1
    chrUn_KI270396v1
    chrUn_KI270388v1
    chrUn_KI270394v1
    chrUn_KI270386v1
    chrUn_KI270391v1
    chrUn_KI270383v1
    chrUn_KI270393v1
    chrUn_KI270384v1
    chrUn_KI270392v1
    chrUn_KI270381v1
    chrUn_KI270385v1
    chrUn_KI270382v1
    chrUn_KI270376v1
    chrUn_KI270374v1
    chrUn_KI270372v1
    chrUn_KI270373v1
    chrUn_KI270375v1
    chrUn_KI270371v1
    chrUn_KI270448v1
    chrUn_KI270521v1
    chrUn_GL000195v1
    chrUn_GL000219v1
    chrUn_GL000220v1
    chrUn_GL000224v1
    chrUn_KI270741v1
    chrUn_GL000226v1
    chrUn_GL000213v1
    chrUn_KI270743v1
    chrUn_KI270744v1
    chrUn_KI270745v1
    chrUn_KI270746v1
    chrUn_KI270747v1
    chrUn_KI270748v1
    chrUn_KI270749v1
    chrUn_KI270750v1
    chrUn_KI270751v1
    chrUn_KI270752v1
    chrUn_KI270753v1
    chrUn_KI270754v1
    chrUn_KI270755v1
    chrUn_KI270756v1
    chrUn_KI270757v1
    chrUn_GL000214v1
    chrUn_KI270742v1
    chrUn_GL000216v2
    chrUn_GL000218v1
    chr1_KI270762v1_alt
    chr1_KI270766v1_alt
    chr1_KI270760v1_alt
    chr1_KI270765v1_alt
    chr1_GL383518v1_alt
    chr1_GL383519v1_alt
    chr1_GL383520v2_alt
    chr1_KI270764v1_alt
    chr1_KI270763v1_alt
    chr1_KI270759v1_alt
    chr1_KI270761v1_alt
    chr2_KI270770v1_alt
    chr2_KI270773v1_alt
    chr2_KI270774v1_alt
    chr2_KI270769v1_alt
    chr2_GL383521v1_alt
    chr2_KI270772v1_alt
    chr2_KI270775v1_alt
    chr2_KI270771v1_alt
    chr2_KI270768v1_alt
    chr2_GL582966v2_alt
    chr2_GL383522v1_alt
    chr2_KI270776v1_alt
    chr2_KI270767v1_alt
    chr3_JH636055v2_alt
    chr3_KI270783v1_alt
    chr3_KI270780v1_alt
    chr3_GL383526v1_alt
    chr3_KI270777v1_alt
    chr3_KI270778v1_alt
    chr3_KI270781v1_alt
    chr3_KI270779v1_alt
    chr3_KI270782v1_alt
    chr3_KI270784v1_alt
    chr4_KI270790v1_alt
    chr4_GL383528v1_alt
    chr4_KI270787v1_alt
    chr4_GL000257v2_alt
    chr4_KI270788v1_alt
    chr4_GL383527v1_alt
    chr4_KI270785v1_alt
    chr4_KI270789v1_alt
    chr4_KI270786v1_alt
    chr5_KI270793v1_alt
    chr5_KI270792v1_alt
    chr5_KI270791v1_alt
    chr5_GL383532v1_alt
    chr5_GL949742v1_alt
    chr5_KI270794v1_alt
    chr5_GL339449v2_alt
    chr5_GL383530v1_alt
    chr5_KI270796v1_alt
    chr5_GL383531v1_alt
    chr5_KI270795v1_alt
    chr6_GL000250v2_alt
    chr6_KI270800v1_alt
    chr6_KI270799v1_alt
    chr6_GL383533v1_alt
    chr6_KI270801v1_alt
    chr6_KI270802v1_alt
    chr6_KB021644v2_alt
    chr6_KI270797v1_alt
    chr6_KI270798v1_alt
    chr7_KI270804v1_alt
    chr7_KI270809v1_alt
    chr7_KI270806v1_alt
    chr7_GL383534v2_alt
    chr7_KI270803v1_alt
    chr7_KI270808v1_alt
    chr7_KI270807v1_alt
    chr7_KI270805v1_alt
    chr8_KI270818v1_alt
    chr8_KI270812v1_alt
    chr8_KI270811v1_alt
    chr8_KI270821v1_alt
    chr8_KI270813v1_alt
    chr8_KI270822v1_alt
    chr8_KI270814v1_alt
    chr8_KI270810v1_alt
    chr8_KI270819v1_alt
    chr8_KI270820v1_alt
    chr8_KI270817v1_alt
    chr8_KI270816v1_alt
    chr8_KI270815v1_alt
    chr9_GL383539v1_alt
    chr9_GL383540v1_alt
    chr9_GL383541v1_alt
    chr9_GL383542v1_alt
    chr9_KI270823v1_alt
    chr10_GL383545v1_alt
    chr10_KI270824v1_alt
    chr10_GL383546v1_alt
    chr10_KI270825v1_alt
    chr11_KI270832v1_alt
    chr11_KI270830v1_alt
    chr11_KI270831v1_alt
    chr11_KI270829v1_alt
    chr11_GL383547v1_alt
    chr11_JH159136v1_alt
    chr11_JH159137v1_alt
    chr11_KI270827v1_alt
    chr11_KI270826v1_alt
    chr12_GL877875v1_alt
    chr12_GL877876v1_alt
    chr12_KI270837v1_alt
    chr12_GL383549v1_alt
    chr12_KI270835v1_alt
    chr12_GL383550v2_alt
    chr12_GL383552v1_alt
    chr12_GL383553v2_alt
    chr12_KI270834v1_alt
    chr12_GL383551v1_alt
    chr12_KI270833v1_alt
    chr12_KI270836v1_alt
    chr13_KI270840v1_alt
    chr13_KI270839v1_alt
    chr13_KI270843v1_alt
    chr13_KI270841v1_alt
    chr13_KI270838v1_alt
    chr13_KI270842v1_alt
    chr14_KI270844v1_alt
    chr14_KI270847v1_alt
    chr14_KI270845v1_alt
    chr14_KI270846v1_alt
    chr15_KI270852v1_alt
    chr15_KI270851v1_alt
    chr15_KI270848v1_alt
    chr15_GL383554v1_alt
    chr15_KI270849v1_alt
    chr15_GL383555v2_alt
    chr15_KI270850v1_alt
    chr16_KI270854v1_alt
    chr16_KI270856v1_alt
    chr16_KI270855v1_alt
    chr16_KI270853v1_alt
    chr16_GL383556v1_alt
    chr16_GL383557v1_alt
    chr17_GL383563v3_alt
    chr17_KI270862v1_alt
    chr17_KI270861v1_alt
    chr17_KI270857v1_alt
    chr17_JH159146v1_alt
    chr17_JH159147v1_alt
    chr17_GL383564v2_alt
    chr17_GL000258v2_alt
    chr17_GL383565v1_alt
    chr17_KI270858v1_alt
    chr17_KI270859v1_alt
    chr17_GL383566v1_alt
    chr17_KI270860v1_alt
    chr18_KI270864v1_alt
    chr18_GL383567v1_alt
    chr18_GL383570v1_alt
    chr18_GL383571v1_alt
    chr18_GL383568v1_alt
    chr18_GL383569v1_alt
    chr18_GL383572v1_alt
    chr18_KI270863v1_alt
    chr19_KI270868v1_alt
    chr19_KI270865v1_alt
    chr19_GL383573v1_alt
    chr19_GL383575v2_alt
    chr19_GL383576v1_alt
    chr19_GL383574v1_alt
    chr19_KI270866v1_alt
    chr19_KI270867v1_alt
    chr19_GL949746v1_alt
    chr20_GL383577v2_alt
    chr20_KI270869v1_alt
    chr20_KI270871v1_alt
    chr20_KI270870v1_alt
    chr21_GL383578v2_alt
    chr21_KI270874v1_alt
    chr21_KI270873v1_alt
    chr21_GL383579v2_alt
    chr21_GL383580v2_alt
    chr21_GL383581v2_alt
    chr21_KI270872v1_alt
    chr22_KI270875v1_alt
    chr22_KI270878v1_alt
    chr22_KI270879v1_alt
    chr22_KI270876v1_alt
    chr22_KI270877v1_alt
    chr22_GL383583v2_alt
    chr22_GL383582v2_alt
    chrX_KI270880v1_alt
    chrX_KI270881v1_alt
    chr19_KI270882v1_alt
    chr19_KI270883v1_alt
    chr19_KI270884v1_alt
    chr19_KI270885v1_alt
    chr19_KI270886v1_alt
    chr19_KI270887v1_alt
    chr19_KI270888v1_alt
    chr19_KI270889v1_alt
    chr19_KI270890v1_alt
    chr19_KI270891v1_alt
    chr1_KI270892v1_alt
    chr2_KI270894v1_alt
    chr2_KI270893v1_alt
    chr3_KI270895v1_alt
    chr4_KI270896v1_alt
    chr5_KI270897v1_alt
    chr5_KI270898v1_alt
    chr6_GL000251v2_alt
    chr7_KI270899v1_alt
    chr8_KI270901v1_alt
    chr8_KI270900v1_alt
    chr11_KI270902v1_alt
    chr11_KI270903v1_alt
    chr12_KI270904v1_alt
    chr15_KI270906v1_alt
    chr15_KI270905v1_alt
    chr17_KI270907v1_alt
    chr17_KI270910v1_alt
    chr17_KI270909v1_alt
    chr17_JH159148v1_alt
    chr17_KI270908v1_alt
    chr18_KI270912v1_alt
    chr18_KI270911v1_alt
    chr19_GL949747v2_alt
    chr22_KB663609v1_alt
    chrX_KI270913v1_alt
    chr19_KI270914v1_alt
    chr19_KI270915v1_alt
    chr19_KI270916v1_alt
    chr19_KI270917v1_alt
    chr19_KI270918v1_alt
    chr19_KI270919v1_alt
    chr19_KI270920v1_alt
    chr19_KI270921v1_alt
    chr19_KI270922v1_alt
    chr19_KI270923v1_alt
    chr3_KI270924v1_alt
    chr4_KI270925v1_alt
    chr6_GL000252v2_alt
    chr8_KI270926v1_alt
    chr11_KI270927v1_alt
    chr19_GL949748v2_alt
    chr22_KI270928v1_alt
    chr19_KI270929v1_alt
    chr19_KI270930v1_alt
    chr19_KI270931v1_alt
    chr19_KI270932v1_alt
    chr19_KI270933v1_alt
    chr19_GL000209v2_alt
    chr3_KI270934v1_alt
    chr6_GL000253v2_alt
    chr19_GL949749v2_alt
    chr3_KI270935v1_alt
    chr6_GL000254v2_alt
    chr19_GL949750v2_alt
    chr3_KI270936v1_alt
    chr6_GL000255v2_alt
    chr19_GL949751v2_alt
    chr3_KI270937v1_alt
    chr6_GL000256v2_alt
    chr19_GL949752v1_alt
    chr6_KI270758v1_alt
    chr19_GL949753v2_alt
    chr19_KI270938v1_alt
    chrEBV
    chrUn_KN707606v1_decoy
    chrUn_KN707607v1_decoy
    chrUn_KN707608v1_decoy
    chrUn_KN707609v1_decoy
    chrUn_KN707610v1_decoy
    chrUn_KN707611v1_decoy
    chrUn_KN707612v1_decoy
    chrUn_KN707613v1_decoy
    chrUn_KN707614v1_decoy
    chrUn_KN707615v1_decoy
    chrUn_KN707616v1_decoy
    chrUn_KN707617v1_decoy
    chrUn_KN707618v1_decoy
    chrUn_KN707619v1_decoy
    chrUn_KN707620v1_decoy
    chrUn_KN707621v1_decoy
    chrUn_KN707622v1_decoy
    chrUn_KN707623v1_decoy
    chrUn_KN707624v1_decoy
    chrUn_KN707625v1_decoy
    chrUn_KN707626v1_decoy
    chrUn_KN707627v1_decoy
    chrUn_KN707628v1_decoy
    chrUn_KN707629v1_decoy
    chrUn_KN707630v1_decoy
    chrUn_KN707631v1_decoy
    chrUn_KN707632v1_decoy
    chrUn_KN707633v1_decoy
    chrUn_KN707634v1_decoy
    chrUn_KN707635v1_decoy
    chrUn_KN707636v1_decoy
    chrUn_KN707637v1_decoy
    chrUn_KN707638v1_decoy
    chrUn_KN707639v1_decoy
    chrUn_KN707640v1_decoy
    chrUn_KN707641v1_decoy
    chrUn_KN707642v1_decoy
    chrUn_KN707643v1_decoy
    chrUn_KN707644v1_decoy
    chrUn_KN707645v1_decoy
    chrUn_KN707646v1_decoy
    chrUn_KN707647v1_decoy
    chrUn_KN707648v1_decoy
    chrUn_KN707649v1_decoy
    chrUn_KN707650v1_decoy
    chrUn_KN707651v1_decoy
    chrUn_KN707652v1_decoy
    chrUn_KN707653v1_decoy
    chrUn_KN707654v1_decoy
    chrUn_KN707655v1_decoy
    chrUn_KN707656v1_decoy
    chrUn_KN707657v1_decoy
    chrUn_KN707658v1_decoy
    chrUn_KN707659v1_decoy
    chrUn_KN707660v1_decoy
    chrUn_KN707661v1_decoy
    chrUn_KN707662v1_decoy
    chrUn_KN707663v1_decoy
    chrUn_KN707664v1_decoy
    chrUn_KN707665v1_decoy
    chrUn_KN707666v1_decoy
    chrUn_KN707667v1_decoy
    chrUn_KN707668v1_decoy
    chrUn_KN707669v1_decoy
    chrUn_KN707670v1_decoy
    chrUn_KN707671v1_decoy
    chrUn_KN707672v1_decoy
    chrUn_KN707673v1_decoy
    chrUn_KN707674v1_decoy
    chrUn_KN707675v1_decoy
    chrUn_KN707676v1_decoy
    chrUn_KN707677v1_decoy
    chrUn_KN707678v1_decoy
    chrUn_KN707679v1_decoy
    chrUn_KN707680v1_decoy
    chrUn_KN707681v1_decoy
    chrUn_KN707682v1_decoy
    chrUn_KN707683v1_decoy
    chrUn_KN707684v1_decoy
    chrUn_KN707685v1_decoy
    chrUn_KN707686v1_decoy
    chrUn_KN707687v1_decoy
    chrUn_KN707688v1_decoy
    chrUn_KN707689v1_decoy
    chrUn_KN707690v1_decoy
    chrUn_KN707691v1_decoy
    chrUn_KN707692v1_decoy
    chrUn_KN707693v1_decoy
    chrUn_KN707694v1_decoy
    chrUn_KN707695v1_decoy
    chrUn_KN707696v1_decoy
    chrUn_KN707697v1_decoy
    chrUn_KN707698v1_decoy
    chrUn_KN707699v1_decoy
    chrUn_KN707700v1_decoy
    chrUn_KN707701v1_decoy
    chrUn_KN707702v1_decoy
    chrUn_KN707703v1_decoy
    chrUn_KN707704v1_decoy
    chrUn_KN707705v1_decoy
    chrUn_KN707706v1_decoy
    chrUn_KN707707v1_decoy
    chrUn_KN707708v1_decoy
    chrUn_KN707709v1_decoy
    chrUn_KN707710v1_decoy
    chrUn_KN707711v1_decoy
    chrUn_KN707712v1_decoy
    chrUn_KN707713v1_decoy
    chrUn_KN707714v1_decoy
    chrUn_KN707715v1_decoy
    chrUn_KN707716v1_decoy
    chrUn_KN707717v1_decoy
    chrUn_KN707718v1_decoy
    chrUn_KN707719v1_decoy
    chrUn_KN707720v1_decoy
    chrUn_KN707721v1_decoy
    chrUn_KN707722v1_decoy
    chrUn_KN707723v1_decoy
    chrUn_KN707724v1_decoy
    chrUn_KN707725v1_decoy
    chrUn_KN707726v1_decoy
    chrUn_KN707727v1_decoy
    chrUn_KN707728v1_decoy
    chrUn_KN707729v1_decoy
    chrUn_KN707730v1_decoy
    chrUn_KN707731v1_decoy
    chrUn_KN707732v1_decoy
    chrUn_KN707733v1_decoy
    chrUn_KN707734v1_decoy
    chrUn_KN707735v1_decoy
    chrUn_KN707736v1_decoy
    chrUn_KN707737v1_decoy
    chrUn_KN707738v1_decoy
    chrUn_KN707739v1_decoy
    chrUn_KN707740v1_decoy
    chrUn_KN707741v1_decoy
    chrUn_KN707742v1_decoy
    chrUn_KN707743v1_decoy
    chrUn_KN707744v1_decoy
    chrUn_KN707745v1_decoy
    chrUn_KN707746v1_decoy
    chrUn_KN707747v1_decoy
    chrUn_KN707748v1_decoy
    chrUn_KN707749v1_decoy
    chrUn_KN707750v1_decoy
    chrUn_KN707751v1_decoy
    chrUn_KN707752v1_decoy
    chrUn_KN707753v1_decoy
    chrUn_KN707754v1_decoy
    chrUn_KN707755v1_decoy
    chrUn_KN707756v1_decoy
    chrUn_KN707757v1_decoy
    chrUn_KN707758v1_decoy
    chrUn_KN707759v1_decoy
    chrUn_KN707760v1_decoy
    chrUn_KN707761v1_decoy
    chrUn_KN707762v1_decoy
    chrUn_KN707763v1_decoy
    chrUn_KN707764v1_decoy
    chrUn_KN707765v1_decoy
    chrUn_KN707766v1_decoy
    chrUn_KN707767v1_decoy
    chrUn_KN707768v1_decoy
    chrUn_KN707769v1_decoy
    chrUn_KN707770v1_decoy
    chrUn_KN707771v1_decoy
    chrUn_KN707772v1_decoy
    chrUn_KN707773v1_decoy
    chrUn_KN707774v1_decoy
    chrUn_KN707775v1_decoy
    chrUn_KN707776v1_decoy
    chrUn_KN707777v1_decoy
    chrUn_KN707778v1_decoy
    chrUn_KN707779v1_decoy
    chrUn_KN707780v1_decoy
    chrUn_KN707781v1_decoy
    chrUn_KN707782v1_decoy
    chrUn_KN707783v1_decoy
    chrUn_KN707784v1_decoy
    chrUn_KN707785v1_decoy
    chrUn_KN707786v1_decoy
    chrUn_KN707787v1_decoy
    chrUn_KN707788v1_decoy
    chrUn_KN707789v1_decoy
    chrUn_KN707790v1_decoy
    chrUn_KN707791v1_decoy
    chrUn_KN707792v1_decoy
    chrUn_KN707793v1_decoy
    chrUn_KN707794v1_decoy
    chrUn_KN707795v1_decoy
    chrUn_KN707796v1_decoy
    chrUn_KN707797v1_decoy
    chrUn_KN707798v1_decoy
    chrUn_KN707799v1_decoy
    chrUn_KN707800v1_decoy
    chrUn_KN707801v1_decoy
    chrUn_KN707802v1_decoy
    chrUn_KN707803v1_decoy
    chrUn_KN707804v1_decoy
    chrUn_KN707805v1_decoy
    chrUn_KN707806v1_decoy
    chrUn_KN707807v1_decoy
    chrUn_KN707808v1_decoy
    chrUn_KN707809v1_decoy
    chrUn_KN707810v1_decoy
    chrUn_KN707811v1_decoy
    chrUn_KN707812v1_decoy
    chrUn_KN707813v1_decoy
    chrUn_KN707814v1_decoy
    chrUn_KN707815v1_decoy
    chrUn_KN707816v1_decoy
    chrUn_KN707817v1_decoy
    chrUn_KN707818v1_decoy
    chrUn_KN707819v1_decoy
    chrUn_KN707820v1_decoy
    chrUn_KN707821v1_decoy
    chrUn_KN707822v1_decoy
    chrUn_KN707823v1_decoy
    chrUn_KN707824v1_decoy
    chrUn_KN707825v1_decoy
    chrUn_KN707826v1_decoy
    chrUn_KN707827v1_decoy
    chrUn_KN707828v1_decoy
    chrUn_KN707829v1_decoy
    chrUn_KN707830v1_decoy
    chrUn_KN707831v1_decoy
    chrUn_KN707832v1_decoy
    chrUn_KN707833v1_decoy
    chrUn_KN707834v1_decoy
    chrUn_KN707835v1_decoy
    chrUn_KN707836v1_decoy
    chrUn_KN707837v1_decoy
    chrUn_KN707838v1_decoy
    chrUn_KN707839v1_decoy
    chrUn_KN707840v1_decoy
    chrUn_KN707841v1_decoy
    chrUn_KN707842v1_decoy
    chrUn_KN707843v1_decoy
    chrUn_KN707844v1_decoy
    chrUn_KN707845v1_decoy
    chrUn_KN707846v1_decoy
    chrUn_KN707847v1_decoy
    chrUn_KN707848v1_decoy
    chrUn_KN707849v1_decoy
    chrUn_KN707850v1_decoy
    chrUn_KN707851v1_decoy
    chrUn_KN707852v1_decoy
    chrUn_KN707853v1_decoy
    chrUn_KN707854v1_decoy
    chrUn_KN707855v1_decoy
    chrUn_KN707856v1_decoy
    chrUn_KN707857v1_decoy
    chrUn_KN707858v1_decoy
    chrUn_KN707859v1_decoy
    chrUn_KN707860v1_decoy
    chrUn_KN707861v1_decoy
    chrUn_KN707862v1_decoy
    chrUn_KN707863v1_decoy
    chrUn_KN707864v1_decoy
    chrUn_KN707865v1_decoy
    chrUn_KN707866v1_decoy
    chrUn_KN707867v1_decoy
    chrUn_KN707868v1_decoy
    chrUn_KN707869v1_decoy
    chrUn_KN707870v1_decoy
    chrUn_KN707871v1_decoy
    chrUn_KN707872v1_decoy
    chrUn_KN707873v1_decoy
    chrUn_KN707874v1_decoy
    chrUn_KN707875v1_decoy
    chrUn_KN707876v1_decoy
    chrUn_KN707877v1_decoy
    chrUn_KN707878v1_decoy
    chrUn_KN707879v1_decoy
    chrUn_KN707880v1_decoy
    chrUn_KN707881v1_decoy
    chrUn_KN707882v1_decoy
    chrUn_KN707883v1_decoy
    chrUn_KN707884v1_decoy
    chrUn_KN707885v1_decoy
    chrUn_KN707886v1_decoy
    chrUn_KN707887v1_decoy
    chrUn_KN707888v1_decoy
    chrUn_KN707889v1_decoy
    chrUn_KN707890v1_decoy
    chrUn_KN707891v1_decoy
    chrUn_KN707892v1_decoy
    chrUn_KN707893v1_decoy
    chrUn_KN707894v1_decoy
    chrUn_KN707895v1_decoy
    chrUn_KN707896v1_decoy
    chrUn_KN707897v1_decoy
    chrUn_KN707898v1_decoy
    chrUn_KN707899v1_decoy
    chrUn_KN707900v1_decoy
    chrUn_KN707901v1_decoy
    chrUn_KN707902v1_decoy
    chrUn_KN707903v1_decoy
    chrUn_KN707904v1_decoy
    chrUn_KN707905v1_decoy
    chrUn_KN707906v1_decoy
    chrUn_KN707907v1_decoy
    chrUn_KN707908v1_decoy
    chrUn_KN707909v1_decoy
    chrUn_KN707910v1_decoy
    chrUn_KN707911v1_decoy
    chrUn_KN707912v1_decoy
    chrUn_KN707913v1_decoy
    chrUn_KN707914v1_decoy
    chrUn_KN707915v1_decoy
    chrUn_KN707916v1_decoy
    chrUn_KN707917v1_decoy
    chrUn_KN707918v1_decoy
    chrUn_KN707919v1_decoy
    chrUn_KN707920v1_decoy
    chrUn_KN707921v1_decoy
    chrUn_KN707922v1_decoy
    chrUn_KN707923v1_decoy
    chrUn_KN707924v1_decoy
    chrUn_KN707925v1_decoy
    chrUn_KN707926v1_decoy
    chrUn_KN707927v1_decoy
    chrUn_KN707928v1_decoy
    chrUn_KN707929v1_decoy
    chrUn_KN707930v1_decoy
    chrUn_KN707931v1_decoy
    chrUn_KN707932v1_decoy
    chrUn_KN707933v1_decoy
    chrUn_KN707934v1_decoy
    chrUn_KN707935v1_decoy
    chrUn_KN707936v1_decoy
    chrUn_KN707937v1_decoy
    chrUn_KN707938v1_decoy
    chrUn_KN707939v1_decoy
    chrUn_KN707940v1_decoy
    chrUn_KN707941v1_decoy
    chrUn_KN707942v1_decoy
    chrUn_KN707943v1_decoy
    chrUn_KN707944v1_decoy
    chrUn_KN707945v1_decoy
    chrUn_KN707946v1_decoy
    chrUn_KN707947v1_decoy
    chrUn_KN707948v1_decoy
    chrUn_KN707949v1_decoy
    chrUn_KN707950v1_decoy
    chrUn_KN707951v1_decoy
    chrUn_KN707952v1_decoy
    chrUn_KN707953v1_decoy
    chrUn_KN707954v1_decoy
    chrUn_KN707955v1_decoy
    chrUn_KN707956v1_decoy
    chrUn_KN707957v1_decoy
    chrUn_KN707958v1_decoy
    chrUn_KN707959v1_decoy
    chrUn_KN707960v1_decoy
    chrUn_KN707961v1_decoy
    chrUn_KN707962v1_decoy
    chrUn_KN707963v1_decoy
    chrUn_KN707964v1_decoy
    chrUn_KN707965v1_decoy
    chrUn_KN707966v1_decoy
    chrUn_KN707967v1_decoy
    chrUn_KN707968v1_decoy
    chrUn_KN707969v1_decoy
    chrUn_KN707970v1_decoy
    chrUn_KN707971v1_decoy
    chrUn_KN707972v1_decoy
    chrUn_KN707973v1_decoy
    chrUn_KN707974v1_decoy
    chrUn_KN707975v1_decoy
    chrUn_KN707976v1_decoy
    chrUn_KN707977v1_decoy
    chrUn_KN707978v1_decoy
    chrUn_KN707979v1_decoy
    chrUn_KN707980v1_decoy
    chrUn_KN707981v1_decoy
    chrUn_KN707982v1_decoy
    chrUn_KN707983v1_decoy
    chrUn_KN707984v1_decoy
    chrUn_KN707985v1_decoy
    chrUn_KN707986v1_decoy
    chrUn_KN707987v1_decoy
    chrUn_KN707988v1_decoy
    chrUn_KN707989v1_decoy
    chrUn_KN707990v1_decoy
    chrUn_KN707991v1_decoy
    chrUn_KN707992v1_decoy
    chrUn_JTFH01000001v1_decoy
    chrUn_JTFH01000002v1_decoy
    chrUn_JTFH01000003v1_decoy
    chrUn_JTFH01000004v1_decoy
    chrUn_JTFH01000005v1_decoy
    chrUn_JTFH01000006v1_decoy
    chrUn_JTFH01000007v1_decoy
    chrUn_JTFH01000008v1_decoy
    chrUn_JTFH01000009v1_decoy
    chrUn_JTFH01000010v1_decoy
    chrUn_JTFH01000011v1_decoy
    chrUn_JTFH01000012v1_decoy
    chrUn_JTFH01000013v1_decoy
    chrUn_JTFH01000014v1_decoy
    chrUn_JTFH01000015v1_decoy
    chrUn_JTFH01000016v1_decoy
    chrUn_JTFH01000017v1_decoy
    chrUn_JTFH01000018v1_decoy
    chrUn_JTFH01000019v1_decoy
    chrUn_JTFH01000020v1_decoy
    chrUn_JTFH01000021v1_decoy
    chrUn_JTFH01000022v1_decoy
    chrUn_JTFH01000023v1_decoy
    chrUn_JTFH01000024v1_decoy
    chrUn_JTFH01000025v1_decoy
    chrUn_JTFH01000026v1_decoy
    chrUn_JTFH01000027v1_decoy
    chrUn_JTFH01000028v1_decoy
    chrUn_JTFH01000029v1_decoy
    chrUn_JTFH01000030v1_decoy
    chrUn_JTFH01000031v1_decoy
    chrUn_JTFH01000032v1_decoy
    chrUn_JTFH01000033v1_decoy
    chrUn_JTFH01000034v1_decoy
    chrUn_JTFH01000035v1_decoy
    chrUn_JTFH01000036v1_decoy
    chrUn_JTFH01000037v1_decoy
    chrUn_JTFH01000038v1_decoy
    chrUn_JTFH01000039v1_decoy
    chrUn_JTFH01000040v1_decoy
    chrUn_JTFH01000041v1_decoy
    chrUn_JTFH01000042v1_decoy
    chrUn_JTFH01000043v1_decoy
    chrUn_JTFH01000044v1_decoy
    chrUn_JTFH01000045v1_decoy
    chrUn_JTFH01000046v1_decoy
    chrUn_JTFH01000047v1_decoy
    chrUn_JTFH01000048v1_decoy
    chrUn_JTFH01000049v1_decoy
    chrUn_JTFH01000050v1_decoy
    chrUn_JTFH01000051v1_decoy
    chrUn_JTFH01000052v1_decoy
    chrUn_JTFH01000053v1_decoy
    chrUn_JTFH01000054v1_decoy
    chrUn_JTFH01000055v1_decoy
    chrUn_JTFH01000056v1_decoy
    chrUn_JTFH01000057v1_decoy
    chrUn_JTFH01000058v1_decoy
    chrUn_JTFH01000059v1_decoy
    chrUn_JTFH01000060v1_decoy
    chrUn_JTFH01000061v1_decoy
    chrUn_JTFH01000062v1_decoy
    chrUn_JTFH01000063v1_decoy
    chrUn_JTFH01000064v1_decoy
    chrUn_JTFH01000065v1_decoy
    chrUn_JTFH01000066v1_decoy
    chrUn_JTFH01000067v1_decoy
    chrUn_JTFH01000068v1_decoy
    chrUn_JTFH01000069v1_decoy
    chrUn_JTFH01000070v1_decoy
    chrUn_JTFH01000071v1_decoy
    chrUn_JTFH01000072v1_decoy
    chrUn_JTFH01000073v1_decoy
    chrUn_JTFH01000074v1_decoy
    chrUn_JTFH01000075v1_decoy
    chrUn_JTFH01000076v1_decoy
    chrUn_JTFH01000077v1_decoy
    chrUn_JTFH01000078v1_decoy
    chrUn_JTFH01000079v1_decoy
    chrUn_JTFH01000080v1_decoy
    chrUn_JTFH01000081v1_decoy
    chrUn_JTFH01000082v1_decoy
    chrUn_JTFH01000083v1_decoy
    chrUn_JTFH01000084v1_decoy
    chrUn_JTFH01000085v1_decoy
    chrUn_JTFH01000086v1_decoy
    chrUn_JTFH01000087v1_decoy
    chrUn_JTFH01000088v1_decoy
    chrUn_JTFH01000089v1_decoy
    chrUn_JTFH01000090v1_decoy
    chrUn_JTFH01000091v1_decoy
    chrUn_JTFH01000092v1_decoy
    chrUn_JTFH01000093v1_decoy
    chrUn_JTFH01000094v1_decoy
    chrUn_JTFH01000095v1_decoy
    chrUn_JTFH01000096v1_decoy
    chrUn_JTFH01000097v1_decoy
    chrUn_JTFH01000098v1_decoy
    chrUn_JTFH01000099v1_decoy
    chrUn_JTFH01000100v1_decoy
    chrUn_JTFH01000101v1_decoy
    chrUn_JTFH01000102v1_decoy
    chrUn_JTFH01000103v1_decoy
    chrUn_JTFH01000104v1_decoy
    chrUn_JTFH01000105v1_decoy
    chrUn_JTFH01000106v1_decoy
    chrUn_JTFH01000107v1_decoy
    chrUn_JTFH01000108v1_decoy
    chrUn_JTFH01000109v1_decoy
    chrUn_JTFH01000110v1_decoy
    chrUn_JTFH01000111v1_decoy
    chrUn_JTFH01000112v1_decoy
    chrUn_JTFH01000113v1_decoy
    chrUn_JTFH01000114v1_decoy
    chrUn_JTFH01000115v1_decoy
    chrUn_JTFH01000116v1_decoy
    chrUn_JTFH01000117v1_decoy
    chrUn_JTFH01000118v1_decoy
    chrUn_JTFH01000119v1_decoy
    chrUn_JTFH01000120v1_decoy
    chrUn_JTFH01000121v1_decoy
    chrUn_JTFH01000122v1_decoy
    chrUn_JTFH01000123v1_decoy
    chrUn_JTFH01000124v1_decoy
    chrUn_JTFH01000125v1_decoy
    chrUn_JTFH01000126v1_decoy
    chrUn_JTFH01000127v1_decoy
    chrUn_JTFH01000128v1_decoy
    chrUn_JTFH01000129v1_decoy
    chrUn_JTFH01000130v1_decoy
    chrUn_JTFH01000131v1_decoy
    chrUn_JTFH01000132v1_decoy
    chrUn_JTFH01000133v1_decoy
    chrUn_JTFH01000134v1_decoy
    chrUn_JTFH01000135v1_decoy
    chrUn_JTFH01000136v1_decoy
    chrUn_JTFH01000137v1_decoy
    chrUn_JTFH01000138v1_decoy
    chrUn_JTFH01000139v1_decoy
    chrUn_JTFH01000140v1_decoy
    chrUn_JTFH01000141v1_decoy
    chrUn_JTFH01000142v1_decoy
    chrUn_JTFH01000143v1_decoy
    chrUn_JTFH01000144v1_decoy
    chrUn_JTFH01000145v1_decoy
    chrUn_JTFH01000146v1_decoy
    chrUn_JTFH01000147v1_decoy
    chrUn_JTFH01000148v1_decoy
    chrUn_JTFH01000149v1_decoy
    chrUn_JTFH01000150v1_decoy
    chrUn_JTFH01000151v1_decoy
    chrUn_JTFH01000152v1_decoy
    chrUn_JTFH01000153v1_decoy
    chrUn_JTFH01000154v1_decoy
    chrUn_JTFH01000155v1_decoy
    chrUn_JTFH01000156v1_decoy
    chrUn_JTFH01000157v1_decoy
    chrUn_JTFH01000158v1_decoy
    chrUn_JTFH01000159v1_decoy
    chrUn_JTFH01000160v1_decoy
    chrUn_JTFH01000161v1_decoy
    chrUn_JTFH01000162v1_decoy
    chrUn_JTFH01000163v1_decoy
    chrUn_JTFH01000164v1_decoy
    chrUn_JTFH01000165v1_decoy
    chrUn_JTFH01000166v1_decoy
    chrUn_JTFH01000167v1_decoy
    chrUn_JTFH01000168v1_decoy
    chrUn_JTFH01000169v1_decoy
    chrUn_JTFH01000170v1_decoy
    chrUn_JTFH01000171v1_decoy
    chrUn_JTFH01000172v1_decoy
    chrUn_JTFH01000173v1_decoy
    chrUn_JTFH01000174v1_decoy
    chrUn_JTFH01000175v1_decoy
    chrUn_JTFH01000176v1_decoy
    chrUn_JTFH01000177v1_decoy
    chrUn_JTFH01000178v1_decoy
    chrUn_JTFH01000179v1_decoy
    chrUn_JTFH01000180v1_decoy
    chrUn_JTFH01000181v1_decoy
    chrUn_JTFH01000182v1_decoy
    chrUn_JTFH01000183v1_decoy
    chrUn_JTFH01000184v1_decoy
    chrUn_JTFH01000185v1_decoy
    chrUn_JTFH01000186v1_decoy
    chrUn_JTFH01000187v1_decoy
    chrUn_JTFH01000188v1_decoy
    chrUn_JTFH01000189v1_decoy
    chrUn_JTFH01000190v1_decoy
    chrUn_JTFH01000191v1_decoy
    chrUn_JTFH01000192v1_decoy
    chrUn_JTFH01000193v1_decoy
    chrUn_JTFH01000194v1_decoy
    chrUn_JTFH01000195v1_decoy
    chrUn_JTFH01000196v1_decoy
    chrUn_JTFH01000197v1_decoy
    chrUn_JTFH01000198v1_decoy
    chrUn_JTFH01000199v1_decoy
    chrUn_JTFH01000200v1_decoy
    chrUn_JTFH01000201v1_decoy
    chrUn_JTFH01000202v1_decoy
    chrUn_JTFH01000203v1_decoy
    chrUn_JTFH01000204v1_decoy
    chrUn_JTFH01000205v1_decoy
    chrUn_JTFH01000206v1_decoy
    chrUn_JTFH01000207v1_decoy
    chrUn_JTFH01000208v1_decoy
    chrUn_JTFH01000209v1_decoy
    chrUn_JTFH01000210v1_decoy
    chrUn_JTFH01000211v1_decoy
    chrUn_JTFH01000212v1_decoy
    chrUn_JTFH01000213v1_decoy
    chrUn_JTFH01000214v1_decoy
    chrUn_JTFH01000215v1_decoy
    chrUn_JTFH01000216v1_decoy
    chrUn_JTFH01000217v1_decoy
    chrUn_JTFH01000218v1_decoy
    chrUn_JTFH01000219v1_decoy
    chrUn_JTFH01000220v1_decoy
    chrUn_JTFH01000221v1_decoy
    chrUn_JTFH01000222v1_decoy
    chrUn_JTFH01000223v1_decoy
    chrUn_JTFH01000224v1_decoy
    chrUn_JTFH01000225v1_decoy
    chrUn_JTFH01000226v1_decoy
    chrUn_JTFH01000227v1_decoy
    chrUn_JTFH01000228v1_decoy
    chrUn_JTFH01000229v1_decoy
    chrUn_JTFH01000230v1_decoy
    chrUn_JTFH01000231v1_decoy
    chrUn_JTFH01000232v1_decoy
    chrUn_JTFH01000233v1_decoy
    chrUn_JTFH01000234v1_decoy
    chrUn_JTFH01000235v1_decoy
    chrUn_JTFH01000236v1_decoy
    chrUn_JTFH01000237v1_decoy
    chrUn_JTFH01000238v1_decoy
    chrUn_JTFH01000239v1_decoy
    chrUn_JTFH01000240v1_decoy
    chrUn_JTFH01000241v1_decoy
    chrUn_JTFH01000242v1_decoy
    chrUn_JTFH01000243v1_decoy
    chrUn_JTFH01000244v1_decoy
    chrUn_JTFH01000245v1_decoy
    chrUn_JTFH01000246v1_decoy
    chrUn_JTFH01000247v1_decoy
    chrUn_JTFH01000248v1_decoy
    chrUn_JTFH01000249v1_decoy
    chrUn_JTFH01000250v1_decoy
    chrUn_JTFH01000251v1_decoy
    chrUn_JTFH01000252v1_decoy
    chrUn_JTFH01000253v1_decoy
    chrUn_JTFH01000254v1_decoy
    chrUn_JTFH01000255v1_decoy
    chrUn_JTFH01000256v1_decoy
    chrUn_JTFH01000257v1_decoy
    chrUn_JTFH01000258v1_decoy
    chrUn_JTFH01000259v1_decoy
    chrUn_JTFH01000260v1_decoy
    chrUn_JTFH01000261v1_decoy
    chrUn_JTFH01000262v1_decoy
    chrUn_JTFH01000263v1_decoy
    chrUn_JTFH01000264v1_decoy
    chrUn_JTFH01000265v1_decoy
    chrUn_JTFH01000266v1_decoy
    chrUn_JTFH01000267v1_decoy
    chrUn_JTFH01000268v1_decoy
    chrUn_JTFH01000269v1_decoy
    chrUn_JTFH01000270v1_decoy
    chrUn_JTFH01000271v1_decoy
    chrUn_JTFH01000272v1_decoy
    chrUn_JTFH01000273v1_decoy
    chrUn_JTFH01000274v1_decoy
    chrUn_JTFH01000275v1_decoy
    chrUn_JTFH01000276v1_decoy
    chrUn_JTFH01000277v1_decoy
    chrUn_JTFH01000278v1_decoy
    chrUn_JTFH01000279v1_decoy
    chrUn_JTFH01000280v1_decoy
    chrUn_JTFH01000281v1_decoy
    chrUn_JTFH01000282v1_decoy
    chrUn_JTFH01000283v1_decoy
    chrUn_JTFH01000284v1_decoy
    chrUn_JTFH01000285v1_decoy
    chrUn_JTFH01000286v1_decoy
    chrUn_JTFH01000287v1_decoy
    chrUn_JTFH01000288v1_decoy
    chrUn_JTFH01000289v1_decoy
    chrUn_JTFH01000290v1_decoy
    chrUn_JTFH01000291v1_decoy
    chrUn_JTFH01000292v1_decoy
    chrUn_JTFH01000293v1_decoy
    chrUn_JTFH01000294v1_decoy
    chrUn_JTFH01000295v1_decoy
    chrUn_JTFH01000296v1_decoy
    chrUn_JTFH01000297v1_decoy
    chrUn_JTFH01000298v1_decoy
    chrUn_JTFH01000299v1_decoy
    chrUn_JTFH01000300v1_decoy
    chrUn_JTFH01000301v1_decoy
    chrUn_JTFH01000302v1_decoy
    chrUn_JTFH01000303v1_decoy
    chrUn_JTFH01000304v1_decoy
    chrUn_JTFH01000305v1_decoy
    chrUn_JTFH01000306v1_decoy
    chrUn_JTFH01000307v1_decoy
    chrUn_JTFH01000308v1_decoy
    chrUn_JTFH01000309v1_decoy
    chrUn_JTFH01000310v1_decoy
    chrUn_JTFH01000311v1_decoy
    chrUn_JTFH01000312v1_decoy
    chrUn_JTFH01000313v1_decoy
    chrUn_JTFH01000314v1_decoy
    chrUn_JTFH01000315v1_decoy
    chrUn_JTFH01000316v1_decoy
    chrUn_JTFH01000317v1_decoy
    chrUn_JTFH01000318v1_decoy
    chrUn_JTFH01000319v1_decoy
    chrUn_JTFH01000320v1_decoy
    chrUn_JTFH01000321v1_decoy
    chrUn_JTFH01000322v1_decoy
    chrUn_JTFH01000323v1_decoy
    chrUn_JTFH01000324v1_decoy
    chrUn_JTFH01000325v1_decoy
    chrUn_JTFH01000326v1_decoy
    chrUn_JTFH01000327v1_decoy
    chrUn_JTFH01000328v1_decoy
    chrUn_JTFH01000329v1_decoy
    chrUn_JTFH01000330v1_decoy
    chrUn_JTFH01000331v1_decoy
    chrUn_JTFH01000332v1_decoy
    chrUn_JTFH01000333v1_decoy
    chrUn_JTFH01000334v1_decoy
    chrUn_JTFH01000335v1_decoy
    chrUn_JTFH01000336v1_decoy
    chrUn_JTFH01000337v1_decoy
    chrUn_JTFH01000338v1_decoy
    chrUn_JTFH01000339v1_decoy
    chrUn_JTFH01000340v1_decoy
    chrUn_JTFH01000341v1_decoy
    chrUn_JTFH01000342v1_decoy
    chrUn_JTFH01000343v1_decoy
    chrUn_JTFH01000344v1_decoy
    chrUn_JTFH01000345v1_decoy
    chrUn_JTFH01000346v1_decoy
    chrUn_JTFH01000347v1_decoy
    chrUn_JTFH01000348v1_decoy
    chrUn_JTFH01000349v1_decoy
    chrUn_JTFH01000350v1_decoy
    chrUn_JTFH01000351v1_decoy
    chrUn_JTFH01000352v1_decoy
    chrUn_JTFH01000353v1_decoy
    chrUn_JTFH01000354v1_decoy
    chrUn_JTFH01000355v1_decoy
    chrUn_JTFH01000356v1_decoy
    chrUn_JTFH01000357v1_decoy
    chrUn_JTFH01000358v1_decoy
    chrUn_JTFH01000359v1_decoy
    chrUn_JTFH01000360v1_decoy
    chrUn_JTFH01000361v1_decoy
    chrUn_JTFH01000362v1_decoy
    chrUn_JTFH01000363v1_decoy
    chrUn_JTFH01000364v1_decoy
    chrUn_JTFH01000365v1_decoy
    chrUn_JTFH01000366v1_decoy
    chrUn_JTFH01000367v1_decoy
    chrUn_JTFH01000368v1_decoy
    chrUn_JTFH01000369v1_decoy
    chrUn_JTFH01000370v1_decoy
    chrUn_JTFH01000371v1_decoy
    chrUn_JTFH01000372v1_decoy
    chrUn_JTFH01000373v1_decoy
    chrUn_JTFH01000374v1_decoy
    chrUn_JTFH01000375v1_decoy
    chrUn_JTFH01000376v1_decoy
    chrUn_JTFH01000377v1_decoy
    chrUn_JTFH01000378v1_decoy
    chrUn_JTFH01000379v1_decoy
    chrUn_JTFH01000380v1_decoy
    chrUn_JTFH01000381v1_decoy
    chrUn_JTFH01000382v1_decoy
    chrUn_JTFH01000383v1_decoy
    chrUn_JTFH01000384v1_decoy
    chrUn_JTFH01000385v1_decoy
    chrUn_JTFH01000386v1_decoy
    chrUn_JTFH01000387v1_decoy
    chrUn_JTFH01000388v1_decoy
    chrUn_JTFH01000389v1_decoy
    chrUn_JTFH01000390v1_decoy
    chrUn_JTFH01000391v1_decoy
    chrUn_JTFH01000392v1_decoy
    chrUn_JTFH01000393v1_decoy
    chrUn_JTFH01000394v1_decoy
    chrUn_JTFH01000395v1_decoy
    chrUn_JTFH01000396v1_decoy
    chrUn_JTFH01000397v1_decoy
    chrUn_JTFH01000398v1_decoy
    chrUn_JTFH01000399v1_decoy
    chrUn_JTFH01000400v1_decoy
    chrUn_JTFH01000401v1_decoy
    chrUn_JTFH01000402v1_decoy
    chrUn_JTFH01000403v1_decoy
    chrUn_JTFH01000404v1_decoy
    chrUn_JTFH01000405v1_decoy
    chrUn_JTFH01000406v1_decoy
    chrUn_JTFH01000407v1_decoy
    chrUn_JTFH01000408v1_decoy
    chrUn_JTFH01000409v1_decoy
    chrUn_JTFH01000410v1_decoy
    chrUn_JTFH01000411v1_decoy
    chrUn_JTFH01000412v1_decoy
    chrUn_JTFH01000413v1_decoy
    chrUn_JTFH01000414v1_decoy
    chrUn_JTFH01000415v1_decoy
    chrUn_JTFH01000416v1_decoy
    chrUn_JTFH01000417v1_decoy
    chrUn_JTFH01000418v1_decoy
    chrUn_JTFH01000419v1_decoy
    chrUn_JTFH01000420v1_decoy
    chrUn_JTFH01000421v1_decoy
    chrUn_JTFH01000422v1_decoy
    chrUn_JTFH01000423v1_decoy
    chrUn_JTFH01000424v1_decoy
    chrUn_JTFH01000425v1_decoy
    chrUn_JTFH01000426v1_decoy
    chrUn_JTFH01000427v1_decoy
    chrUn_JTFH01000428v1_decoy
    chrUn_JTFH01000429v1_decoy
    chrUn_JTFH01000430v1_decoy
    chrUn_JTFH01000431v1_decoy
    chrUn_JTFH01000432v1_decoy
    chrUn_JTFH01000433v1_decoy
    chrUn_JTFH01000434v1_decoy
    chrUn_JTFH01000435v1_decoy
    chrUn_JTFH01000436v1_decoy
    chrUn_JTFH01000437v1_decoy
    chrUn_JTFH01000438v1_decoy
    chrUn_JTFH01000439v1_decoy
    chrUn_JTFH01000440v1_decoy
    chrUn_JTFH01000441v1_decoy
    chrUn_JTFH01000442v1_decoy
    chrUn_JTFH01000443v1_decoy
    chrUn_JTFH01000444v1_decoy
    chrUn_JTFH01000445v1_decoy
    chrUn_JTFH01000446v1_decoy
    chrUn_JTFH01000447v1_decoy
    chrUn_JTFH01000448v1_decoy
    chrUn_JTFH01000449v1_decoy
    chrUn_JTFH01000450v1_decoy
    chrUn_JTFH01000451v1_decoy
    chrUn_JTFH01000452v1_decoy
    chrUn_JTFH01000453v1_decoy
    chrUn_JTFH01000454v1_decoy
    chrUn_JTFH01000455v1_decoy
    chrUn_JTFH01000456v1_decoy
    chrUn_JTFH01000457v1_decoy
    chrUn_JTFH01000458v1_decoy
    chrUn_JTFH01000459v1_decoy
    chrUn_JTFH01000460v1_decoy
    chrUn_JTFH01000461v1_decoy
    chrUn_JTFH01000462v1_decoy
    chrUn_JTFH01000463v1_decoy
    chrUn_JTFH01000464v1_decoy
    chrUn_JTFH01000465v1_decoy
    chrUn_JTFH01000466v1_decoy
    chrUn_JTFH01000467v1_decoy
    chrUn_JTFH01000468v1_decoy
    chrUn_JTFH01000469v1_decoy
    chrUn_JTFH01000470v1_decoy
    chrUn_JTFH01000471v1_decoy
    chrUn_JTFH01000472v1_decoy
    chrUn_JTFH01000473v1_decoy
    chrUn_JTFH01000474v1_decoy
    chrUn_JTFH01000475v1_decoy
    chrUn_JTFH01000476v1_decoy
    chrUn_JTFH01000477v1_decoy
    chrUn_JTFH01000478v1_decoy
    chrUn_JTFH01000479v1_decoy
    chrUn_JTFH01000480v1_decoy
    chrUn_JTFH01000481v1_decoy
    chrUn_JTFH01000482v1_decoy
    chrUn_JTFH01000483v1_decoy
    chrUn_JTFH01000484v1_decoy
    chrUn_JTFH01000485v1_decoy
    chrUn_JTFH01000486v1_decoy
    chrUn_JTFH01000487v1_decoy
    chrUn_JTFH01000488v1_decoy
    chrUn_JTFH01000489v1_decoy
    chrUn_JTFH01000490v1_decoy
    chrUn_JTFH01000491v1_decoy
    chrUn_JTFH01000492v1_decoy
    chrUn_JTFH01000493v1_decoy
    chrUn_JTFH01000494v1_decoy
    chrUn_JTFH01000495v1_decoy
    chrUn_JTFH01000496v1_decoy
    chrUn_JTFH01000497v1_decoy
    chrUn_JTFH01000498v1_decoy
    chrUn_JTFH01000499v1_decoy
    chrUn_JTFH01000500v1_decoy
    chrUn_JTFH01000501v1_decoy
    chrUn_JTFH01000502v1_decoy
    chrUn_JTFH01000503v1_decoy
    chrUn_JTFH01000504v1_decoy
    chrUn_JTFH01000505v1_decoy
    chrUn_JTFH01000506v1_decoy
    chrUn_JTFH01000507v1_decoy
    chrUn_JTFH01000508v1_decoy
    chrUn_JTFH01000509v1_decoy
    chrUn_JTFH01000510v1_decoy
    chrUn_JTFH01000511v1_decoy
    chrUn_JTFH01000512v1_decoy
    chrUn_JTFH01000513v1_decoy
    chrUn_JTFH01000514v1_decoy
    chrUn_JTFH01000515v1_decoy
    chrUn_JTFH01000516v1_decoy
    chrUn_JTFH01000517v1_decoy
    chrUn_JTFH01000518v1_decoy
    chrUn_JTFH01000519v1_decoy
    chrUn_JTFH01000520v1_decoy
    chrUn_JTFH01000521v1_decoy
    chrUn_JTFH01000522v1_decoy
    chrUn_JTFH01000523v1_decoy
    chrUn_JTFH01000524v1_decoy
    chrUn_JTFH01000525v1_decoy
    chrUn_JTFH01000526v1_decoy
    chrUn_JTFH01000527v1_decoy
    chrUn_JTFH01000528v1_decoy
    chrUn_JTFH01000529v1_decoy
    chrUn_JTFH01000530v1_decoy
    chrUn_JTFH01000531v1_decoy
    chrUn_JTFH01000532v1_decoy
    chrUn_JTFH01000533v1_decoy
    chrUn_JTFH01000534v1_decoy
    chrUn_JTFH01000535v1_decoy
    chrUn_JTFH01000536v1_decoy
    chrUn_JTFH01000537v1_decoy
    chrUn_JTFH01000538v1_decoy
    chrUn_JTFH01000539v1_decoy
    chrUn_JTFH01000540v1_decoy
    chrUn_JTFH01000541v1_decoy
    chrUn_JTFH01000542v1_decoy
    chrUn_JTFH01000543v1_decoy
    chrUn_JTFH01000544v1_decoy
    chrUn_JTFH01000545v1_decoy
    chrUn_JTFH01000546v1_decoy
    chrUn_JTFH01000547v1_decoy
    chrUn_JTFH01000548v1_decoy
    chrUn_JTFH01000549v1_decoy
    chrUn_JTFH01000550v1_decoy
    chrUn_JTFH01000551v1_decoy
    chrUn_JTFH01000552v1_decoy
    chrUn_JTFH01000553v1_decoy
    chrUn_JTFH01000554v1_decoy
    chrUn_JTFH01000555v1_decoy
    chrUn_JTFH01000556v1_decoy
    chrUn_JTFH01000557v1_decoy
    chrUn_JTFH01000558v1_decoy
    chrUn_JTFH01000559v1_decoy
    chrUn_JTFH01000560v1_decoy
    chrUn_JTFH01000561v1_decoy
    chrUn_JTFH01000562v1_decoy
    chrUn_JTFH01000563v1_decoy
    chrUn_JTFH01000564v1_decoy
    chrUn_JTFH01000565v1_decoy
    chrUn_JTFH01000566v1_decoy
    chrUn_JTFH01000567v1_decoy
    chrUn_JTFH01000568v1_decoy
    chrUn_JTFH01000569v1_decoy
    chrUn_JTFH01000570v1_decoy
    chrUn_JTFH01000571v1_decoy
    chrUn_JTFH01000572v1_decoy
    chrUn_JTFH01000573v1_decoy
    chrUn_JTFH01000574v1_decoy
    chrUn_JTFH01000575v1_decoy
    chrUn_JTFH01000576v1_decoy
    chrUn_JTFH01000577v1_decoy
    chrUn_JTFH01000578v1_decoy
    chrUn_JTFH01000579v1_decoy
    chrUn_JTFH01000580v1_decoy
    chrUn_JTFH01000581v1_decoy
    chrUn_JTFH01000582v1_decoy
    chrUn_JTFH01000583v1_decoy
    chrUn_JTFH01000584v1_decoy
    chrUn_JTFH01000585v1_decoy
    chrUn_JTFH01000586v1_decoy
    chrUn_JTFH01000587v1_decoy
    chrUn_JTFH01000588v1_decoy
    chrUn_JTFH01000589v1_decoy
    chrUn_JTFH01000590v1_decoy
    chrUn_JTFH01000591v1_decoy
    chrUn_JTFH01000592v1_decoy
    chrUn_JTFH01000593v1_decoy
    chrUn_JTFH01000594v1_decoy
    chrUn_JTFH01000595v1_decoy
    chrUn_JTFH01000596v1_decoy
    chrUn_JTFH01000597v1_decoy
    chrUn_JTFH01000598v1_decoy
    chrUn_JTFH01000599v1_decoy
    chrUn_JTFH01000600v1_decoy
    chrUn_JTFH01000601v1_decoy
    chrUn_JTFH01000602v1_decoy
    chrUn_JTFH01000603v1_decoy
    chrUn_JTFH01000604v1_decoy
    chrUn_JTFH01000605v1_decoy
    chrUn_JTFH01000606v1_decoy
    chrUn_JTFH01000607v1_decoy
    chrUn_JTFH01000608v1_decoy
    chrUn_JTFH01000609v1_decoy
    chrUn_JTFH01000610v1_decoy
    chrUn_JTFH01000611v1_decoy
    chrUn_JTFH01000612v1_decoy
    chrUn_JTFH01000613v1_decoy
    chrUn_JTFH01000614v1_decoy
    chrUn_JTFH01000615v1_decoy
    chrUn_JTFH01000616v1_decoy
    chrUn_JTFH01000617v1_decoy
    chrUn_JTFH01000618v1_decoy
    chrUn_JTFH01000619v1_decoy
    chrUn_JTFH01000620v1_decoy
    chrUn_JTFH01000621v1_decoy
    chrUn_JTFH01000622v1_decoy
    chrUn_JTFH01000623v1_decoy
    chrUn_JTFH01000624v1_decoy
    chrUn_JTFH01000625v1_decoy
    chrUn_JTFH01000626v1_decoy
    chrUn_JTFH01000627v1_decoy
    chrUn_JTFH01000628v1_decoy
    chrUn_JTFH01000629v1_decoy
    chrUn_JTFH01000630v1_decoy
    chrUn_JTFH01000631v1_decoy
    chrUn_JTFH01000632v1_decoy
    chrUn_JTFH01000633v1_decoy
    chrUn_JTFH01000634v1_decoy
    chrUn_JTFH01000635v1_decoy
    chrUn_JTFH01000636v1_decoy
    chrUn_JTFH01000637v1_decoy
    chrUn_JTFH01000638v1_decoy
    chrUn_JTFH01000639v1_decoy
    chrUn_JTFH01000640v1_decoy
    chrUn_JTFH01000641v1_decoy
    chrUn_JTFH01000642v1_decoy
    chrUn_JTFH01000643v1_decoy
    chrUn_JTFH01000644v1_decoy
    chrUn_JTFH01000645v1_decoy
    chrUn_JTFH01000646v1_decoy
    chrUn_JTFH01000647v1_decoy
    chrUn_JTFH01000648v1_decoy
    chrUn_JTFH01000649v1_decoy
    chrUn_JTFH01000650v1_decoy
    chrUn_JTFH01000651v1_decoy
    chrUn_JTFH01000652v1_decoy
    chrUn_JTFH01000653v1_decoy
    chrUn_JTFH01000654v1_decoy
    chrUn_JTFH01000655v1_decoy
    chrUn_JTFH01000656v1_decoy
    chrUn_JTFH01000657v1_decoy
    chrUn_JTFH01000658v1_decoy
    chrUn_JTFH01000659v1_decoy
    chrUn_JTFH01000660v1_decoy
    chrUn_JTFH01000661v1_decoy
    chrUn_JTFH01000662v1_decoy
    chrUn_JTFH01000663v1_decoy
    chrUn_JTFH01000664v1_decoy
    chrUn_JTFH01000665v1_decoy
    chrUn_JTFH01000666v1_decoy
    chrUn_JTFH01000667v1_decoy
    chrUn_JTFH01000668v1_decoy
    chrUn_JTFH01000669v1_decoy
    chrUn_JTFH01000670v1_decoy
    chrUn_JTFH01000671v1_decoy
    chrUn_JTFH01000672v1_decoy
    chrUn_JTFH01000673v1_decoy
    chrUn_JTFH01000674v1_decoy
    chrUn_JTFH01000675v1_decoy
    chrUn_JTFH01000676v1_decoy
    chrUn_JTFH01000677v1_decoy
    chrUn_JTFH01000678v1_decoy
    chrUn_JTFH01000679v1_decoy
    chrUn_JTFH01000680v1_decoy
    chrUn_JTFH01000681v1_decoy
    chrUn_JTFH01000682v1_decoy
    chrUn_JTFH01000683v1_decoy
    chrUn_JTFH01000684v1_decoy
    chrUn_JTFH01000685v1_decoy
    chrUn_JTFH01000686v1_decoy
    chrUn_JTFH01000687v1_decoy
    chrUn_JTFH01000688v1_decoy
    chrUn_JTFH01000689v1_decoy
    chrUn_JTFH01000690v1_decoy
    chrUn_JTFH01000691v1_decoy
    chrUn_JTFH01000692v1_decoy
    chrUn_JTFH01000693v1_decoy
    chrUn_JTFH01000694v1_decoy
    chrUn_JTFH01000695v1_decoy
    chrUn_JTFH01000696v1_decoy
    chrUn_JTFH01000697v1_decoy
    chrUn_JTFH01000698v1_decoy
    chrUn_JTFH01000699v1_decoy
    chrUn_JTFH01000700v1_decoy
    chrUn_JTFH01000701v1_decoy
    chrUn_JTFH01000702v1_decoy
    chrUn_JTFH01000703v1_decoy
    chrUn_JTFH01000704v1_decoy
    chrUn_JTFH01000705v1_decoy
    chrUn_JTFH01000706v1_decoy
    chrUn_JTFH01000707v1_decoy
    chrUn_JTFH01000708v1_decoy
    chrUn_JTFH01000709v1_decoy
    chrUn_JTFH01000710v1_decoy
    chrUn_JTFH01000711v1_decoy
    chrUn_JTFH01000712v1_decoy
    chrUn_JTFH01000713v1_decoy
    chrUn_JTFH01000714v1_decoy
    chrUn_JTFH01000715v1_decoy
    chrUn_JTFH01000716v1_decoy
    chrUn_JTFH01000717v1_decoy
    chrUn_JTFH01000718v1_decoy
    chrUn_JTFH01000719v1_decoy
    chrUn_JTFH01000720v1_decoy
    chrUn_JTFH01000721v1_decoy
    chrUn_JTFH01000722v1_decoy
    chrUn_JTFH01000723v1_decoy
    chrUn_JTFH01000724v1_decoy
    chrUn_JTFH01000725v1_decoy
    chrUn_JTFH01000726v1_decoy
    chrUn_JTFH01000727v1_decoy
    chrUn_JTFH01000728v1_decoy
    chrUn_JTFH01000729v1_decoy
    chrUn_JTFH01000730v1_decoy
    chrUn_JTFH01000731v1_decoy
    chrUn_JTFH01000732v1_decoy
    chrUn_JTFH01000733v1_decoy
    chrUn_JTFH01000734v1_decoy
    chrUn_JTFH01000735v1_decoy
    chrUn_JTFH01000736v1_decoy
    chrUn_JTFH01000737v1_decoy
    chrUn_JTFH01000738v1_decoy
    chrUn_JTFH01000739v1_decoy
    chrUn_JTFH01000740v1_decoy
    chrUn_JTFH01000741v1_decoy
    chrUn_JTFH01000742v1_decoy
    chrUn_JTFH01000743v1_decoy
    chrUn_JTFH01000744v1_decoy
    chrUn_JTFH01000745v1_decoy
    chrUn_JTFH01000746v1_decoy
    chrUn_JTFH01000747v1_decoy
    chrUn_JTFH01000748v1_decoy
    chrUn_JTFH01000749v1_decoy
    chrUn_JTFH01000750v1_decoy
    chrUn_JTFH01000751v1_decoy
    chrUn_JTFH01000752v1_decoy
    chrUn_JTFH01000753v1_decoy
    chrUn_JTFH01000754v1_decoy
    chrUn_JTFH01000755v1_decoy
    chrUn_JTFH01000756v1_decoy
    chrUn_JTFH01000757v1_decoy
    chrUn_JTFH01000758v1_decoy
    chrUn_JTFH01000759v1_decoy
    chrUn_JTFH01000760v1_decoy
    chrUn_JTFH01000761v1_decoy
    chrUn_JTFH01000762v1_decoy
    chrUn_JTFH01000763v1_decoy
    chrUn_JTFH01000764v1_decoy
    chrUn_JTFH01000765v1_decoy
    chrUn_JTFH01000766v1_decoy
    chrUn_JTFH01000767v1_decoy
    chrUn_JTFH01000768v1_decoy
    chrUn_JTFH01000769v1_decoy
    chrUn_JTFH01000770v1_decoy
    chrUn_JTFH01000771v1_decoy
    chrUn_JTFH01000772v1_decoy
    chrUn_JTFH01000773v1_decoy
    chrUn_JTFH01000774v1_decoy
    chrUn_JTFH01000775v1_decoy
    chrUn_JTFH01000776v1_decoy
    chrUn_JTFH01000777v1_decoy
    chrUn_JTFH01000778v1_decoy
    chrUn_JTFH01000779v1_decoy
    chrUn_JTFH01000780v1_decoy
    chrUn_JTFH01000781v1_decoy
    chrUn_JTFH01000782v1_decoy
    chrUn_JTFH01000783v1_decoy
    chrUn_JTFH01000784v1_decoy
    chrUn_JTFH01000785v1_decoy
    chrUn_JTFH01000786v1_decoy
    chrUn_JTFH01000787v1_decoy
    chrUn_JTFH01000788v1_decoy
    chrUn_JTFH01000789v1_decoy
    chrUn_JTFH01000790v1_decoy
    chrUn_JTFH01000791v1_decoy
    chrUn_JTFH01000792v1_decoy
    chrUn_JTFH01000793v1_decoy
    chrUn_JTFH01000794v1_decoy
    chrUn_JTFH01000795v1_decoy
    chrUn_JTFH01000796v1_decoy
    chrUn_JTFH01000797v1_decoy
    chrUn_JTFH01000798v1_decoy
    chrUn_JTFH01000799v1_decoy
    chrUn_JTFH01000800v1_decoy
    chrUn_JTFH01000801v1_decoy
    chrUn_JTFH01000802v1_decoy
    chrUn_JTFH01000803v1_decoy
    chrUn_JTFH01000804v1_decoy
    chrUn_JTFH01000805v1_decoy
    chrUn_JTFH01000806v1_decoy
    chrUn_JTFH01000807v1_decoy
    chrUn_JTFH01000808v1_decoy
    chrUn_JTFH01000809v1_decoy
    chrUn_JTFH01000810v1_decoy
    chrUn_JTFH01000811v1_decoy
    chrUn_JTFH01000812v1_decoy
    chrUn_JTFH01000813v1_decoy
    chrUn_JTFH01000814v1_decoy
    chrUn_JTFH01000815v1_decoy
    chrUn_JTFH01000816v1_decoy
    chrUn_JTFH01000817v1_decoy
    chrUn_JTFH01000818v1_decoy
    chrUn_JTFH01000819v1_decoy
    chrUn_JTFH01000820v1_decoy
    chrUn_JTFH01000821v1_decoy
    chrUn_JTFH01000822v1_decoy
    chrUn_JTFH01000823v1_decoy
    chrUn_JTFH01000824v1_decoy
    chrUn_JTFH01000825v1_decoy
    chrUn_JTFH01000826v1_decoy
    chrUn_JTFH01000827v1_decoy
    chrUn_JTFH01000828v1_decoy
    chrUn_JTFH01000829v1_decoy
    chrUn_JTFH01000830v1_decoy
    chrUn_JTFH01000831v1_decoy
    chrUn_JTFH01000832v1_decoy
    chrUn_JTFH01000833v1_decoy
    chrUn_JTFH01000834v1_decoy
    chrUn_JTFH01000835v1_decoy
    chrUn_JTFH01000836v1_decoy
    chrUn_JTFH01000837v1_decoy
    chrUn_JTFH01000838v1_decoy
    chrUn_JTFH01000839v1_decoy
    chrUn_JTFH01000840v1_decoy
    chrUn_JTFH01000841v1_decoy
    chrUn_JTFH01000842v1_decoy
    chrUn_JTFH01000843v1_decoy
    chrUn_JTFH01000844v1_decoy
    chrUn_JTFH01000845v1_decoy
    chrUn_JTFH01000846v1_decoy
    chrUn_JTFH01000847v1_decoy
    chrUn_JTFH01000848v1_decoy
    chrUn_JTFH01000849v1_decoy
    chrUn_JTFH01000850v1_decoy
    chrUn_JTFH01000851v1_decoy
    chrUn_JTFH01000852v1_decoy
    chrUn_JTFH01000853v1_decoy
    chrUn_JTFH01000854v1_decoy
    chrUn_JTFH01000855v1_decoy
    chrUn_JTFH01000856v1_decoy
    chrUn_JTFH01000857v1_decoy
    chrUn_JTFH01000858v1_decoy
    chrUn_JTFH01000859v1_decoy
    chrUn_JTFH01000860v1_decoy
    chrUn_JTFH01000861v1_decoy
    chrUn_JTFH01000862v1_decoy
    chrUn_JTFH01000863v1_decoy
    chrUn_JTFH01000864v1_decoy
    chrUn_JTFH01000865v1_decoy
    chrUn_JTFH01000866v1_decoy
    chrUn_JTFH01000867v1_decoy
    chrUn_JTFH01000868v1_decoy
    chrUn_JTFH01000869v1_decoy
    chrUn_JTFH01000870v1_decoy
    chrUn_JTFH01000871v1_decoy
    chrUn_JTFH01000872v1_decoy
    chrUn_JTFH01000873v1_decoy
    chrUn_JTFH01000874v1_decoy
    chrUn_JTFH01000875v1_decoy
    chrUn_JTFH01000876v1_decoy
    chrUn_JTFH01000877v1_decoy
    chrUn_JTFH01000878v1_decoy
    chrUn_JTFH01000879v1_decoy
    chrUn_JTFH01000880v1_decoy
    chrUn_JTFH01000881v1_decoy
    chrUn_JTFH01000882v1_decoy
    chrUn_JTFH01000883v1_decoy
    chrUn_JTFH01000884v1_decoy
    chrUn_JTFH01000885v1_decoy
    chrUn_JTFH01000886v1_decoy
    chrUn_JTFH01000887v1_decoy
    chrUn_JTFH01000888v1_decoy
    chrUn_JTFH01000889v1_decoy
    chrUn_JTFH01000890v1_decoy
    chrUn_JTFH01000891v1_decoy
    chrUn_JTFH01000892v1_decoy
    chrUn_JTFH01000893v1_decoy
    chrUn_JTFH01000894v1_decoy
    chrUn_JTFH01000895v1_decoy
    chrUn_JTFH01000896v1_decoy
    chrUn_JTFH01000897v1_decoy
    chrUn_JTFH01000898v1_decoy
    chrUn_JTFH01000899v1_decoy
    chrUn_JTFH01000900v1_decoy
    chrUn_JTFH01000901v1_decoy
    chrUn_JTFH01000902v1_decoy
    chrUn_JTFH01000903v1_decoy
    chrUn_JTFH01000904v1_decoy
    chrUn_JTFH01000905v1_decoy
    chrUn_JTFH01000906v1_decoy
    chrUn_JTFH01000907v1_decoy
    chrUn_JTFH01000908v1_decoy
    chrUn_JTFH01000909v1_decoy
    chrUn_JTFH01000910v1_decoy
    chrUn_JTFH01000911v1_decoy
    chrUn_JTFH01000912v1_decoy
    chrUn_JTFH01000913v1_decoy
    chrUn_JTFH01000914v1_decoy
    chrUn_JTFH01000915v1_decoy
    chrUn_JTFH01000916v1_decoy
    chrUn_JTFH01000917v1_decoy
    chrUn_JTFH01000918v1_decoy
    chrUn_JTFH01000919v1_decoy
    chrUn_JTFH01000920v1_decoy
    chrUn_JTFH01000921v1_decoy
    chrUn_JTFH01000922v1_decoy
    chrUn_JTFH01000923v1_decoy
    chrUn_JTFH01000924v1_decoy
    chrUn_JTFH01000925v1_decoy
    chrUn_JTFH01000926v1_decoy
    chrUn_JTFH01000927v1_decoy
    chrUn_JTFH01000928v1_decoy
    chrUn_JTFH01000929v1_decoy
    chrUn_JTFH01000930v1_decoy
    chrUn_JTFH01000931v1_decoy
    chrUn_JTFH01000932v1_decoy
    chrUn_JTFH01000933v1_decoy
    chrUn_JTFH01000934v1_decoy
    chrUn_JTFH01000935v1_decoy
    chrUn_JTFH01000936v1_decoy
    chrUn_JTFH01000937v1_decoy
    chrUn_JTFH01000938v1_decoy
    chrUn_JTFH01000939v1_decoy
    chrUn_JTFH01000940v1_decoy
    chrUn_JTFH01000941v1_decoy
    chrUn_JTFH01000942v1_decoy
    chrUn_JTFH01000943v1_decoy
    chrUn_JTFH01000944v1_decoy
    chrUn_JTFH01000945v1_decoy
    chrUn_JTFH01000946v1_decoy
    chrUn_JTFH01000947v1_decoy
    chrUn_JTFH01000948v1_decoy
    chrUn_JTFH01000949v1_decoy
    chrUn_JTFH01000950v1_decoy
    chrUn_JTFH01000951v1_decoy
    chrUn_JTFH01000952v1_decoy
    chrUn_JTFH01000953v1_decoy
    chrUn_JTFH01000954v1_decoy
    chrUn_JTFH01000955v1_decoy
    chrUn_JTFH01000956v1_decoy
    chrUn_JTFH01000957v1_decoy
    chrUn_JTFH01000958v1_decoy
    chrUn_JTFH01000959v1_decoy
    chrUn_JTFH01000960v1_decoy
    chrUn_JTFH01000961v1_decoy
    chrUn_JTFH01000962v1_decoy
    chrUn_JTFH01000963v1_decoy
    chrUn_JTFH01000964v1_decoy
    chrUn_JTFH01000965v1_decoy
    chrUn_JTFH01000966v1_decoy
    chrUn_JTFH01000967v1_decoy
    chrUn_JTFH01000968v1_decoy
    chrUn_JTFH01000969v1_decoy
    chrUn_JTFH01000970v1_decoy
    chrUn_JTFH01000971v1_decoy
    chrUn_JTFH01000972v1_decoy
    chrUn_JTFH01000973v1_decoy
    chrUn_JTFH01000974v1_decoy
    chrUn_JTFH01000975v1_decoy
    chrUn_JTFH01000976v1_decoy
    chrUn_JTFH01000977v1_decoy
    chrUn_JTFH01000978v1_decoy
    chrUn_JTFH01000979v1_decoy
    chrUn_JTFH01000980v1_decoy
    chrUn_JTFH01000981v1_decoy
    chrUn_JTFH01000982v1_decoy
    chrUn_JTFH01000983v1_decoy
    chrUn_JTFH01000984v1_decoy
    chrUn_JTFH01000985v1_decoy
    chrUn_JTFH01000986v1_decoy
    chrUn_JTFH01000987v1_decoy
    chrUn_JTFH01000988v1_decoy
    chrUn_JTFH01000989v1_decoy
    chrUn_JTFH01000990v1_decoy
    chrUn_JTFH01000991v1_decoy
    chrUn_JTFH01000992v1_decoy
    chrUn_JTFH01000993v1_decoy
    chrUn_JTFH01000994v1_decoy
    chrUn_JTFH01000995v1_decoy
    chrUn_JTFH01000996v1_decoy
    chrUn_JTFH01000997v1_decoy
    chrUn_JTFH01000998v1_decoy
    chrUn_JTFH01000999v1_decoy
    chrUn_JTFH01001000v1_decoy
    chrUn_JTFH01001001v1_decoy
    chrUn_JTFH01001002v1_decoy
    chrUn_JTFH01001003v1_decoy
    chrUn_JTFH01001004v1_decoy
    chrUn_JTFH01001005v1_decoy
    chrUn_JTFH01001006v1_decoy
    chrUn_JTFH01001007v1_decoy
    chrUn_JTFH01001008v1_decoy
    chrUn_JTFH01001009v1_decoy
    chrUn_JTFH01001010v1_decoy
    chrUn_JTFH01001011v1_decoy
    chrUn_JTFH01001012v1_decoy
    chrUn_JTFH01001013v1_decoy
    chrUn_JTFH01001014v1_decoy
    chrUn_JTFH01001015v1_decoy
    chrUn_JTFH01001016v1_decoy
    chrUn_JTFH01001017v1_decoy
    chrUn_JTFH01001018v1_decoy
    chrUn_JTFH01001019v1_decoy
    chrUn_JTFH01001020v1_decoy
    chrUn_JTFH01001021v1_decoy
    chrUn_JTFH01001022v1_decoy
    chrUn_JTFH01001023v1_decoy
    chrUn_JTFH01001024v1_decoy
    chrUn_JTFH01001025v1_decoy
    chrUn_JTFH01001026v1_decoy
    chrUn_JTFH01001027v1_decoy
    chrUn_JTFH01001028v1_decoy
    chrUn_JTFH01001029v1_decoy
    chrUn_JTFH01001030v1_decoy
    chrUn_JTFH01001031v1_decoy
    chrUn_JTFH01001032v1_decoy
    chrUn_JTFH01001033v1_decoy
    chrUn_JTFH01001034v1_decoy
    chrUn_JTFH01001035v1_decoy
    chrUn_JTFH01001036v1_decoy
    chrUn_JTFH01001037v1_decoy
    chrUn_JTFH01001038v1_decoy
    chrUn_JTFH01001039v1_decoy
    chrUn_JTFH01001040v1_decoy
    chrUn_JTFH01001041v1_decoy
    chrUn_JTFH01001042v1_decoy
    chrUn_JTFH01001043v1_decoy
    chrUn_JTFH01001044v1_decoy
    chrUn_JTFH01001045v1_decoy
    chrUn_JTFH01001046v1_decoy
    chrUn_JTFH01001047v1_decoy
    chrUn_JTFH01001048v1_decoy
    chrUn_JTFH01001049v1_decoy
    chrUn_JTFH01001050v1_decoy
    chrUn_JTFH01001051v1_decoy
    chrUn_JTFH01001052v1_decoy
    chrUn_JTFH01001053v1_decoy
    chrUn_JTFH01001054v1_decoy
    chrUn_JTFH01001055v1_decoy
    chrUn_JTFH01001056v1_decoy
    chrUn_JTFH01001057v1_decoy
    chrUn_JTFH01001058v1_decoy
    chrUn_JTFH01001059v1_decoy
    chrUn_JTFH01001060v1_decoy
    chrUn_JTFH01001061v1_decoy
    chrUn_JTFH01001062v1_decoy
    chrUn_JTFH01001063v1_decoy
    chrUn_JTFH01001064v1_decoy
    chrUn_JTFH01001065v1_decoy
    chrUn_JTFH01001066v1_decoy
    chrUn_JTFH01001067v1_decoy
    chrUn_JTFH01001068v1_decoy
    chrUn_JTFH01001069v1_decoy
    chrUn_JTFH01001070v1_decoy
    chrUn_JTFH01001071v1_decoy
    chrUn_JTFH01001072v1_decoy
    chrUn_JTFH01001073v1_decoy
    chrUn_JTFH01001074v1_decoy
    chrUn_JTFH01001075v1_decoy
    chrUn_JTFH01001076v1_decoy
    chrUn_JTFH01001077v1_decoy
    chrUn_JTFH01001078v1_decoy
    chrUn_JTFH01001079v1_decoy
    chrUn_JTFH01001080v1_decoy
    chrUn_JTFH01001081v1_decoy
    chrUn_JTFH01001082v1_decoy
    chrUn_JTFH01001083v1_decoy
    chrUn_JTFH01001084v1_decoy
    chrUn_JTFH01001085v1_decoy
    chrUn_JTFH01001086v1_decoy
    chrUn_JTFH01001087v1_decoy
    chrUn_JTFH01001088v1_decoy
    chrUn_JTFH01001089v1_decoy
    chrUn_JTFH01001090v1_decoy
    chrUn_JTFH01001091v1_decoy
    chrUn_JTFH01001092v1_decoy
    chrUn_JTFH01001093v1_decoy
    chrUn_JTFH01001094v1_decoy
    chrUn_JTFH01001095v1_decoy
    chrUn_JTFH01001096v1_decoy
    chrUn_JTFH01001097v1_decoy
    chrUn_JTFH01001098v1_decoy
    chrUn_JTFH01001099v1_decoy
    chrUn_JTFH01001100v1_decoy
    chrUn_JTFH01001101v1_decoy
    chrUn_JTFH01001102v1_decoy
    chrUn_JTFH01001103v1_decoy
    chrUn_JTFH01001104v1_decoy
    chrUn_JTFH01001105v1_decoy
    chrUn_JTFH01001106v1_decoy
    chrUn_JTFH01001107v1_decoy
    chrUn_JTFH01001108v1_decoy
    chrUn_JTFH01001109v1_decoy
    chrUn_JTFH01001110v1_decoy
    chrUn_JTFH01001111v1_decoy
    chrUn_JTFH01001112v1_decoy
    chrUn_JTFH01001113v1_decoy
    chrUn_JTFH01001114v1_decoy
    chrUn_JTFH01001115v1_decoy
    chrUn_JTFH01001116v1_decoy
    chrUn_JTFH01001117v1_decoy
    chrUn_JTFH01001118v1_decoy
    chrUn_JTFH01001119v1_decoy
    chrUn_JTFH01001120v1_decoy
    chrUn_JTFH01001121v1_decoy
    chrUn_JTFH01001122v1_decoy
    chrUn_JTFH01001123v1_decoy
    chrUn_JTFH01001124v1_decoy
    chrUn_JTFH01001125v1_decoy
    chrUn_JTFH01001126v1_decoy
    chrUn_JTFH01001127v1_decoy
    chrUn_JTFH01001128v1_decoy
    chrUn_JTFH01001129v1_decoy
    chrUn_JTFH01001130v1_decoy
    chrUn_JTFH01001131v1_decoy
    chrUn_JTFH01001132v1_decoy
    chrUn_JTFH01001133v1_decoy
    chrUn_JTFH01001134v1_decoy
    chrUn_JTFH01001135v1_decoy
    chrUn_JTFH01001136v1_decoy
    chrUn_JTFH01001137v1_decoy
    chrUn_JTFH01001138v1_decoy
    chrUn_JTFH01001139v1_decoy
    chrUn_JTFH01001140v1_decoy
    chrUn_JTFH01001141v1_decoy
    chrUn_JTFH01001142v1_decoy
    chrUn_JTFH01001143v1_decoy
    chrUn_JTFH01001144v1_decoy
    chrUn_JTFH01001145v1_decoy
    chrUn_JTFH01001146v1_decoy
    chrUn_JTFH01001147v1_decoy
    chrUn_JTFH01001148v1_decoy
    chrUn_JTFH01001149v1_decoy
    chrUn_JTFH01001150v1_decoy
    chrUn_JTFH01001151v1_decoy
    chrUn_JTFH01001152v1_decoy
    chrUn_JTFH01001153v1_decoy
    chrUn_JTFH01001154v1_decoy
    chrUn_JTFH01001155v1_decoy
    chrUn_JTFH01001156v1_decoy
    chrUn_JTFH01001157v1_decoy
    chrUn_JTFH01001158v1_decoy
    chrUn_JTFH01001159v1_decoy
    chrUn_JTFH01001160v1_decoy
    chrUn_JTFH01001161v1_decoy
    chrUn_JTFH01001162v1_decoy
    chrUn_JTFH01001163v1_decoy
    chrUn_JTFH01001164v1_decoy
    chrUn_JTFH01001165v1_decoy
    chrUn_JTFH01001166v1_decoy
    chrUn_JTFH01001167v1_decoy
    chrUn_JTFH01001168v1_decoy
    chrUn_JTFH01001169v1_decoy
    chrUn_JTFH01001170v1_decoy
    chrUn_JTFH01001171v1_decoy
    chrUn_JTFH01001172v1_decoy
    chrUn_JTFH01001173v1_decoy
    chrUn_JTFH01001174v1_decoy
    chrUn_JTFH01001175v1_decoy
    chrUn_JTFH01001176v1_decoy
    chrUn_JTFH01001177v1_decoy
    chrUn_JTFH01001178v1_decoy
    chrUn_JTFH01001179v1_decoy
    chrUn_JTFH01001180v1_decoy
    chrUn_JTFH01001181v1_decoy
    chrUn_JTFH01001182v1_decoy
    chrUn_JTFH01001183v1_decoy
    chrUn_JTFH01001184v1_decoy
    chrUn_JTFH01001185v1_decoy
    chrUn_JTFH01001186v1_decoy
    chrUn_JTFH01001187v1_decoy
    chrUn_JTFH01001188v1_decoy
    chrUn_JTFH01001189v1_decoy
    chrUn_JTFH01001190v1_decoy
    chrUn_JTFH01001191v1_decoy
    chrUn_JTFH01001192v1_decoy
    chrUn_JTFH01001193v1_decoy
    chrUn_JTFH01001194v1_decoy
    chrUn_JTFH01001195v1_decoy
    chrUn_JTFH01001196v1_decoy
    chrUn_JTFH01001197v1_decoy
    chrUn_JTFH01001198v1_decoy
    chrUn_JTFH01001199v1_decoy
    chrUn_JTFH01001200v1_decoy
    chrUn_JTFH01001201v1_decoy
    chrUn_JTFH01001202v1_decoy
    chrUn_JTFH01001203v1_decoy
    chrUn_JTFH01001204v1_decoy
    chrUn_JTFH01001205v1_decoy
    chrUn_JTFH01001206v1_decoy
    chrUn_JTFH01001207v1_decoy
    chrUn_JTFH01001208v1_decoy
    chrUn_JTFH01001209v1_decoy
    chrUn_JTFH01001210v1_decoy
    chrUn_JTFH01001211v1_decoy
    chrUn_JTFH01001212v1_decoy
    chrUn_JTFH01001213v1_decoy
    chrUn_JTFH01001214v1_decoy
    chrUn_JTFH01001215v1_decoy
    chrUn_JTFH01001216v1_decoy
    chrUn_JTFH01001217v1_decoy
    chrUn_JTFH01001218v1_decoy
    chrUn_JTFH01001219v1_decoy
    chrUn_JTFH01001220v1_decoy
    chrUn_JTFH01001221v1_decoy
    chrUn_JTFH01001222v1_decoy
    chrUn_JTFH01001223v1_decoy
    chrUn_JTFH01001224v1_decoy
    chrUn_JTFH01001225v1_decoy
    chrUn_JTFH01001226v1_decoy
    chrUn_JTFH01001227v1_decoy
    chrUn_JTFH01001228v1_decoy
    chrUn_JTFH01001229v1_decoy
    chrUn_JTFH01001230v1_decoy
    chrUn_JTFH01001231v1_decoy
    chrUn_JTFH01001232v1_decoy
    chrUn_JTFH01001233v1_decoy
    chrUn_JTFH01001234v1_decoy
    chrUn_JTFH01001235v1_decoy
    chrUn_JTFH01001236v1_decoy
    chrUn_JTFH01001237v1_decoy
    chrUn_JTFH01001238v1_decoy
    chrUn_JTFH01001239v1_decoy
    chrUn_JTFH01001240v1_decoy
    chrUn_JTFH01001241v1_decoy
    chrUn_JTFH01001242v1_decoy
    chrUn_JTFH01001243v1_decoy
    chrUn_JTFH01001244v1_decoy
    chrUn_JTFH01001245v1_decoy
    chrUn_JTFH01001246v1_decoy
    chrUn_JTFH01001247v1_decoy
    chrUn_JTFH01001248v1_decoy
    chrUn_JTFH01001249v1_decoy
    chrUn_JTFH01001250v1_decoy
    chrUn_JTFH01001251v1_decoy
    chrUn_JTFH01001252v1_decoy
    chrUn_JTFH01001253v1_decoy
    chrUn_JTFH01001254v1_decoy
    chrUn_JTFH01001255v1_decoy
    chrUn_JTFH01001256v1_decoy
    chrUn_JTFH01001257v1_decoy
    chrUn_JTFH01001258v1_decoy
    chrUn_JTFH01001259v1_decoy
    chrUn_JTFH01001260v1_decoy
    chrUn_JTFH01001261v1_decoy
    chrUn_JTFH01001262v1_decoy
    chrUn_JTFH01001263v1_decoy
    chrUn_JTFH01001264v1_decoy
    chrUn_JTFH01001265v1_decoy
    chrUn_JTFH01001266v1_decoy
    chrUn_JTFH01001267v1_decoy
    chrUn_JTFH01001268v1_decoy
    chrUn_JTFH01001269v1_decoy
    chrUn_JTFH01001270v1_decoy
    chrUn_JTFH01001271v1_decoy
    chrUn_JTFH01001272v1_decoy
    chrUn_JTFH01001273v1_decoy
    chrUn_JTFH01001274v1_decoy
    chrUn_JTFH01001275v1_decoy
    chrUn_JTFH01001276v1_decoy
    chrUn_JTFH01001277v1_decoy
    chrUn_JTFH01001278v1_decoy
    chrUn_JTFH01001279v1_decoy
    chrUn_JTFH01001280v1_decoy
    chrUn_JTFH01001281v1_decoy
    chrUn_JTFH01001282v1_decoy
    chrUn_JTFH01001283v1_decoy
    chrUn_JTFH01001284v1_decoy
    chrUn_JTFH01001285v1_decoy
    chrUn_JTFH01001286v1_decoy
    chrUn_JTFH01001287v1_decoy
    chrUn_JTFH01001288v1_decoy
    chrUn_JTFH01001289v1_decoy
    chrUn_JTFH01001290v1_decoy
    chrUn_JTFH01001291v1_decoy
    chrUn_JTFH01001292v1_decoy
    chrUn_JTFH01001293v1_decoy
    chrUn_JTFH01001294v1_decoy
    chrUn_JTFH01001295v1_decoy
    chrUn_JTFH01001296v1_decoy
    chrUn_JTFH01001297v1_decoy
    chrUn_JTFH01001298v1_decoy
    chrUn_JTFH01001299v1_decoy
    chrUn_JTFH01001300v1_decoy
    chrUn_JTFH01001301v1_decoy
    chrUn_JTFH01001302v1_decoy
    chrUn_JTFH01001303v1_decoy
    chrUn_JTFH01001304v1_decoy
    chrUn_JTFH01001305v1_decoy
    chrUn_JTFH01001306v1_decoy
    chrUn_JTFH01001307v1_decoy
    chrUn_JTFH01001308v1_decoy
    chrUn_JTFH01001309v1_decoy
    chrUn_JTFH01001310v1_decoy
    chrUn_JTFH01001311v1_decoy
    chrUn_JTFH01001312v1_decoy
    chrUn_JTFH01001313v1_decoy
    chrUn_JTFH01001314v1_decoy
    chrUn_JTFH01001315v1_decoy
    chrUn_JTFH01001316v1_decoy
    chrUn_JTFH01001317v1_decoy
    chrUn_JTFH01001318v1_decoy
    chrUn_JTFH01001319v1_decoy
    chrUn_JTFH01001320v1_decoy
    chrUn_JTFH01001321v1_decoy
    chrUn_JTFH01001322v1_decoy
    chrUn_JTFH01001323v1_decoy
    chrUn_JTFH01001324v1_decoy
    chrUn_JTFH01001325v1_decoy
    chrUn_JTFH01001326v1_decoy
    chrUn_JTFH01001327v1_decoy
    chrUn_JTFH01001328v1_decoy
    chrUn_JTFH01001329v1_decoy
    chrUn_JTFH01001330v1_decoy
    chrUn_JTFH01001331v1_decoy
    chrUn_JTFH01001332v1_decoy
    chrUn_JTFH01001333v1_decoy
    chrUn_JTFH01001334v1_decoy
    chrUn_JTFH01001335v1_decoy
    chrUn_JTFH01001336v1_decoy
    chrUn_JTFH01001337v1_decoy
    chrUn_JTFH01001338v1_decoy
    chrUn_JTFH01001339v1_decoy
    chrUn_JTFH01001340v1_decoy
    chrUn_JTFH01001341v1_decoy
    chrUn_JTFH01001342v1_decoy
    chrUn_JTFH01001343v1_decoy
    chrUn_JTFH01001344v1_decoy
    chrUn_JTFH01001345v1_decoy
    chrUn_JTFH01001346v1_decoy
    chrUn_JTFH01001347v1_decoy
    chrUn_JTFH01001348v1_decoy
    chrUn_JTFH01001349v1_decoy
    chrUn_JTFH01001350v1_decoy
    chrUn_JTFH01001351v1_decoy
    chrUn_JTFH01001352v1_decoy
    chrUn_JTFH01001353v1_decoy
    chrUn_JTFH01001354v1_decoy
    chrUn_JTFH01001355v1_decoy
    chrUn_JTFH01001356v1_decoy
    chrUn_JTFH01001357v1_decoy
    chrUn_JTFH01001358v1_decoy
    chrUn_JTFH01001359v1_decoy
    chrUn_JTFH01001360v1_decoy
    chrUn_JTFH01001361v1_decoy
    chrUn_JTFH01001362v1_decoy
    chrUn_JTFH01001363v1_decoy
    chrUn_JTFH01001364v1_decoy
    chrUn_JTFH01001365v1_decoy
    chrUn_JTFH01001366v1_decoy
    chrUn_JTFH01001367v1_decoy
    chrUn_JTFH01001368v1_decoy
    chrUn_JTFH01001369v1_decoy
    chrUn_JTFH01001370v1_decoy
    chrUn_JTFH01001371v1_decoy
    chrUn_JTFH01001372v1_decoy
    chrUn_JTFH01001373v1_decoy
    chrUn_JTFH01001374v1_decoy
    chrUn_JTFH01001375v1_decoy
    chrUn_JTFH01001376v1_decoy
    chrUn_JTFH01001377v1_decoy
    chrUn_JTFH01001378v1_decoy
    chrUn_JTFH01001379v1_decoy
    chrUn_JTFH01001380v1_decoy
    chrUn_JTFH01001381v1_decoy
    chrUn_JTFH01001382v1_decoy
    chrUn_JTFH01001383v1_decoy
    chrUn_JTFH01001384v1_decoy
    chrUn_JTFH01001385v1_decoy
    chrUn_JTFH01001386v1_decoy
    chrUn_JTFH01001387v1_decoy
    chrUn_JTFH01001388v1_decoy
    chrUn_JTFH01001389v1_decoy
    chrUn_JTFH01001390v1_decoy
    chrUn_JTFH01001391v1_decoy
    chrUn_JTFH01001392v1_decoy
    chrUn_JTFH01001393v1_decoy
    chrUn_JTFH01001394v1_decoy
    chrUn_JTFH01001395v1_decoy
    chrUn_JTFH01001396v1_decoy
    chrUn_JTFH01001397v1_decoy
    chrUn_JTFH01001398v1_decoy
    chrUn_JTFH01001399v1_decoy
    chrUn_JTFH01001400v1_decoy
    chrUn_JTFH01001401v1_decoy
    chrUn_JTFH01001402v1_decoy
    chrUn_JTFH01001403v1_decoy
    chrUn_JTFH01001404v1_decoy
    chrUn_JTFH01001405v1_decoy
    chrUn_JTFH01001406v1_decoy
    chrUn_JTFH01001407v1_decoy
    chrUn_JTFH01001408v1_decoy
    chrUn_JTFH01001409v1_decoy
    chrUn_JTFH01001410v1_decoy
    chrUn_JTFH01001411v1_decoy
    chrUn_JTFH01001412v1_decoy
    chrUn_JTFH01001413v1_decoy
    chrUn_JTFH01001414v1_decoy
    chrUn_JTFH01001415v1_decoy
    chrUn_JTFH01001416v1_decoy
    chrUn_JTFH01001417v1_decoy
    chrUn_JTFH01001418v1_decoy
    chrUn_JTFH01001419v1_decoy
    chrUn_JTFH01001420v1_decoy
    chrUn_JTFH01001421v1_decoy
    chrUn_JTFH01001422v1_decoy
    chrUn_JTFH01001423v1_decoy
    chrUn_JTFH01001424v1_decoy
    chrUn_JTFH01001425v1_decoy
    chrUn_JTFH01001426v1_decoy
    chrUn_JTFH01001427v1_decoy
    chrUn_JTFH01001428v1_decoy
    chrUn_JTFH01001429v1_decoy
    chrUn_JTFH01001430v1_decoy
    chrUn_JTFH01001431v1_decoy
    chrUn_JTFH01001432v1_decoy
    chrUn_JTFH01001433v1_decoy
    chrUn_JTFH01001434v1_decoy
    chrUn_JTFH01001435v1_decoy
    chrUn_JTFH01001436v1_decoy
    chrUn_JTFH01001437v1_decoy
    chrUn_JTFH01001438v1_decoy
    chrUn_JTFH01001439v1_decoy
    chrUn_JTFH01001440v1_decoy
    chrUn_JTFH01001441v1_decoy
    chrUn_JTFH01001442v1_decoy
    chrUn_JTFH01001443v1_decoy
    chrUn_JTFH01001444v1_decoy
    chrUn_JTFH01001445v1_decoy
    chrUn_JTFH01001446v1_decoy
    chrUn_JTFH01001447v1_decoy
    chrUn_JTFH01001448v1_decoy
    chrUn_JTFH01001449v1_decoy
    chrUn_JTFH01001450v1_decoy
    chrUn_JTFH01001451v1_decoy
    chrUn_JTFH01001452v1_decoy
    chrUn_JTFH01001453v1_decoy
    chrUn_JTFH01001454v1_decoy
    chrUn_JTFH01001455v1_decoy
    chrUn_JTFH01001456v1_decoy
    chrUn_JTFH01001457v1_decoy
    chrUn_JTFH01001458v1_decoy
    chrUn_JTFH01001459v1_decoy
    chrUn_JTFH01001460v1_decoy
    chrUn_JTFH01001461v1_decoy
    chrUn_JTFH01001462v1_decoy
    chrUn_JTFH01001463v1_decoy
    chrUn_JTFH01001464v1_decoy
    chrUn_JTFH01001465v1_decoy
    chrUn_JTFH01001466v1_decoy
    chrUn_JTFH01001467v1_decoy
    chrUn_JTFH01001468v1_decoy
    chrUn_JTFH01001469v1_decoy
    chrUn_JTFH01001470v1_decoy
    chrUn_JTFH01001471v1_decoy
    chrUn_JTFH01001472v1_decoy
    chrUn_JTFH01001473v1_decoy
    chrUn_JTFH01001474v1_decoy
    chrUn_JTFH01001475v1_decoy
    chrUn_JTFH01001476v1_decoy
    chrUn_JTFH01001477v1_decoy
    chrUn_JTFH01001478v1_decoy
    chrUn_JTFH01001479v1_decoy
    chrUn_JTFH01001480v1_decoy
    chrUn_JTFH01001481v1_decoy
    chrUn_JTFH01001482v1_decoy
    chrUn_JTFH01001483v1_decoy
    chrUn_JTFH01001484v1_decoy
    chrUn_JTFH01001485v1_decoy
    chrUn_JTFH01001486v1_decoy
    chrUn_JTFH01001487v1_decoy
    chrUn_JTFH01001488v1_decoy
    chrUn_JTFH01001489v1_decoy
    chrUn_JTFH01001490v1_decoy
    chrUn_JTFH01001491v1_decoy
    chrUn_JTFH01001492v1_decoy
    chrUn_JTFH01001493v1_decoy
    chrUn_JTFH01001494v1_decoy
    chrUn_JTFH01001495v1_decoy
    chrUn_JTFH01001496v1_decoy
    chrUn_JTFH01001497v1_decoy
    chrUn_JTFH01001498v1_decoy
    chrUn_JTFH01001499v1_decoy
    chrUn_JTFH01001500v1_decoy
    chrUn_JTFH01001501v1_decoy
    chrUn_JTFH01001502v1_decoy
    chrUn_JTFH01001503v1_decoy
    chrUn_JTFH01001504v1_decoy
    chrUn_JTFH01001505v1_decoy
    chrUn_JTFH01001506v1_decoy
    chrUn_JTFH01001507v1_decoy
    chrUn_JTFH01001508v1_decoy
    chrUn_JTFH01001509v1_decoy
    chrUn_JTFH01001510v1_decoy
    chrUn_JTFH01001511v1_decoy
    chrUn_JTFH01001512v1_decoy
    chrUn_JTFH01001513v1_decoy
    chrUn_JTFH01001514v1_decoy
    chrUn_JTFH01001515v1_decoy
    chrUn_JTFH01001516v1_decoy
    chrUn_JTFH01001517v1_decoy
    chrUn_JTFH01001518v1_decoy
    chrUn_JTFH01001519v1_decoy
    chrUn_JTFH01001520v1_decoy
    chrUn_JTFH01001521v1_decoy
    chrUn_JTFH01001522v1_decoy
    chrUn_JTFH01001523v1_decoy
    chrUn_JTFH01001524v1_decoy
    chrUn_JTFH01001525v1_decoy
    chrUn_JTFH01001526v1_decoy
    chrUn_JTFH01001527v1_decoy
    chrUn_JTFH01001528v1_decoy
    chrUn_JTFH01001529v1_decoy
    chrUn_JTFH01001530v1_decoy
    chrUn_JTFH01001531v1_decoy
    chrUn_JTFH01001532v1_decoy
    chrUn_JTFH01001533v1_decoy
    chrUn_JTFH01001534v1_decoy
    chrUn_JTFH01001535v1_decoy
    chrUn_JTFH01001536v1_decoy
    chrUn_JTFH01001537v1_decoy
    chrUn_JTFH01001538v1_decoy
    chrUn_JTFH01001539v1_decoy
    chrUn_JTFH01001540v1_decoy
    chrUn_JTFH01001541v1_decoy
    chrUn_JTFH01001542v1_decoy
    chrUn_JTFH01001543v1_decoy
    chrUn_JTFH01001544v1_decoy
    chrUn_JTFH01001545v1_decoy
    chrUn_JTFH01001546v1_decoy
    chrUn_JTFH01001547v1_decoy
    chrUn_JTFH01001548v1_decoy
    chrUn_JTFH01001549v1_decoy
    chrUn_JTFH01001550v1_decoy
    chrUn_JTFH01001551v1_decoy
    chrUn_JTFH01001552v1_decoy
    chrUn_JTFH01001553v1_decoy
    chrUn_JTFH01001554v1_decoy
    chrUn_JTFH01001555v1_decoy
    chrUn_JTFH01001556v1_decoy
    chrUn_JTFH01001557v1_decoy
    chrUn_JTFH01001558v1_decoy
    chrUn_JTFH01001559v1_decoy
    chrUn_JTFH01001560v1_decoy
    chrUn_JTFH01001561v1_decoy
    chrUn_JTFH01001562v1_decoy
    chrUn_JTFH01001563v1_decoy
    chrUn_JTFH01001564v1_decoy
    chrUn_JTFH01001565v1_decoy
    chrUn_JTFH01001566v1_decoy
    chrUn_JTFH01001567v1_decoy
    chrUn_JTFH01001568v1_decoy
    chrUn_JTFH01001569v1_decoy
    chrUn_JTFH01001570v1_decoy
    chrUn_JTFH01001571v1_decoy
    chrUn_JTFH01001572v1_decoy
    chrUn_JTFH01001573v1_decoy
    chrUn_JTFH01001574v1_decoy
    chrUn_JTFH01001575v1_decoy
    chrUn_JTFH01001576v1_decoy
    chrUn_JTFH01001577v1_decoy
    chrUn_JTFH01001578v1_decoy
    chrUn_JTFH01001579v1_decoy
    chrUn_JTFH01001580v1_decoy
    chrUn_JTFH01001581v1_decoy
    chrUn_JTFH01001582v1_decoy
    chrUn_JTFH01001583v1_decoy
    chrUn_JTFH01001584v1_decoy
    chrUn_JTFH01001585v1_decoy
    chrUn_JTFH01001586v1_decoy
    chrUn_JTFH01001587v1_decoy
    chrUn_JTFH01001588v1_decoy
    chrUn_JTFH01001589v1_decoy
    chrUn_JTFH01001590v1_decoy
    chrUn_JTFH01001591v1_decoy
    chrUn_JTFH01001592v1_decoy
    chrUn_JTFH01001593v1_decoy
    chrUn_JTFH01001594v1_decoy
    chrUn_JTFH01001595v1_decoy
    chrUn_JTFH01001596v1_decoy
    chrUn_JTFH01001597v1_decoy
    chrUn_JTFH01001598v1_decoy
    chrUn_JTFH01001599v1_decoy
    chrUn_JTFH01001600v1_decoy
    chrUn_JTFH01001601v1_decoy
    chrUn_JTFH01001602v1_decoy
    chrUn_JTFH01001603v1_decoy
    chrUn_JTFH01001604v1_decoy
    chrUn_JTFH01001605v1_decoy
    chrUn_JTFH01001606v1_decoy
    chrUn_JTFH01001607v1_decoy
    chrUn_JTFH01001608v1_decoy
    chrUn_JTFH01001609v1_decoy
    chrUn_JTFH01001610v1_decoy
    chrUn_JTFH01001611v1_decoy
    chrUn_JTFH01001612v1_decoy
    chrUn_JTFH01001613v1_decoy
    chrUn_JTFH01001614v1_decoy
    chrUn_JTFH01001615v1_decoy
    chrUn_JTFH01001616v1_decoy
    chrUn_JTFH01001617v1_decoy
    chrUn_JTFH01001618v1_decoy
    chrUn_JTFH01001619v1_decoy
    chrUn_JTFH01001620v1_decoy
    chrUn_JTFH01001621v1_decoy
    chrUn_JTFH01001622v1_decoy
    chrUn_JTFH01001623v1_decoy
    chrUn_JTFH01001624v1_decoy
    chrUn_JTFH01001625v1_decoy
    chrUn_JTFH01001626v1_decoy
    chrUn_JTFH01001627v1_decoy
    chrUn_JTFH01001628v1_decoy
    chrUn_JTFH01001629v1_decoy
    chrUn_JTFH01001630v1_decoy
    chrUn_JTFH01001631v1_decoy
    chrUn_JTFH01001632v1_decoy
    chrUn_JTFH01001633v1_decoy
    chrUn_JTFH01001634v1_decoy
    chrUn_JTFH01001635v1_decoy
    chrUn_JTFH01001636v1_decoy
    chrUn_JTFH01001637v1_decoy
    chrUn_JTFH01001638v1_decoy
    chrUn_JTFH01001639v1_decoy
    chrUn_JTFH01001640v1_decoy
    chrUn_JTFH01001641v1_decoy
    chrUn_JTFH01001642v1_decoy
    chrUn_JTFH01001643v1_decoy
    chrUn_JTFH01001644v1_decoy
    chrUn_JTFH01001645v1_decoy
    chrUn_JTFH01001646v1_decoy
    chrUn_JTFH01001647v1_decoy
    chrUn_JTFH01001648v1_decoy
    chrUn_JTFH01001649v1_decoy
    chrUn_JTFH01001650v1_decoy
    chrUn_JTFH01001651v1_decoy
    chrUn_JTFH01001652v1_decoy
    chrUn_JTFH01001653v1_decoy
    chrUn_JTFH01001654v1_decoy
    chrUn_JTFH01001655v1_decoy
    chrUn_JTFH01001656v1_decoy
    chrUn_JTFH01001657v1_decoy
    chrUn_JTFH01001658v1_decoy
    chrUn_JTFH01001659v1_decoy
    chrUn_JTFH01001660v1_decoy
    chrUn_JTFH01001661v1_decoy
    chrUn_JTFH01001662v1_decoy
    chrUn_JTFH01001663v1_decoy
    chrUn_JTFH01001664v1_decoy
    chrUn_JTFH01001665v1_decoy
    chrUn_JTFH01001666v1_decoy
    chrUn_JTFH01001667v1_decoy
    chrUn_JTFH01001668v1_decoy
    chrUn_JTFH01001669v1_decoy
    chrUn_JTFH01001670v1_decoy
    chrUn_JTFH01001671v1_decoy
    chrUn_JTFH01001672v1_decoy
    chrUn_JTFH01001673v1_decoy
    chrUn_JTFH01001674v1_decoy
    chrUn_JTFH01001675v1_decoy
    chrUn_JTFH01001676v1_decoy
    chrUn_JTFH01001677v1_decoy
    chrUn_JTFH01001678v1_decoy
    chrUn_JTFH01001679v1_decoy
    chrUn_JTFH01001680v1_decoy
    chrUn_JTFH01001681v1_decoy
    chrUn_JTFH01001682v1_decoy
    chrUn_JTFH01001683v1_decoy
    chrUn_JTFH01001684v1_decoy
    chrUn_JTFH01001685v1_decoy
    chrUn_JTFH01001686v1_decoy
    chrUn_JTFH01001687v1_decoy
    chrUn_JTFH01001688v1_decoy
    chrUn_JTFH01001689v1_decoy
    chrUn_JTFH01001690v1_decoy
    chrUn_JTFH01001691v1_decoy
    chrUn_JTFH01001692v1_decoy
    chrUn_JTFH01001693v1_decoy
    chrUn_JTFH01001694v1_decoy
    chrUn_JTFH01001695v1_decoy
    chrUn_JTFH01001696v1_decoy
    chrUn_JTFH01001697v1_decoy
    chrUn_JTFH01001698v1_decoy
    chrUn_JTFH01001699v1_decoy
    chrUn_JTFH01001700v1_decoy
    chrUn_JTFH01001701v1_decoy
    chrUn_JTFH01001702v1_decoy
    chrUn_JTFH01001703v1_decoy
    chrUn_JTFH01001704v1_decoy
    chrUn_JTFH01001705v1_decoy
    chrUn_JTFH01001706v1_decoy
    chrUn_JTFH01001707v1_decoy
    chrUn_JTFH01001708v1_decoy
    chrUn_JTFH01001709v1_decoy
    chrUn_JTFH01001710v1_decoy
    chrUn_JTFH01001711v1_decoy
    chrUn_JTFH01001712v1_decoy
    chrUn_JTFH01001713v1_decoy
    chrUn_JTFH01001714v1_decoy
    chrUn_JTFH01001715v1_decoy
    chrUn_JTFH01001716v1_decoy
    chrUn_JTFH01001717v1_decoy
    chrUn_JTFH01001718v1_decoy
    chrUn_JTFH01001719v1_decoy
    chrUn_JTFH01001720v1_decoy
    chrUn_JTFH01001721v1_decoy
    chrUn_JTFH01001722v1_decoy
    chrUn_JTFH01001723v1_decoy
    chrUn_JTFH01001724v1_decoy
    chrUn_JTFH01001725v1_decoy
    chrUn_JTFH01001726v1_decoy
    chrUn_JTFH01001727v1_decoy
    chrUn_JTFH01001728v1_decoy
    chrUn_JTFH01001729v1_decoy
    chrUn_JTFH01001730v1_decoy
    chrUn_JTFH01001731v1_decoy
    chrUn_JTFH01001732v1_decoy
    chrUn_JTFH01001733v1_decoy
    chrUn_JTFH01001734v1_decoy
    chrUn_JTFH01001735v1_decoy
    chrUn_JTFH01001736v1_decoy
    chrUn_JTFH01001737v1_decoy
    chrUn_JTFH01001738v1_decoy
    chrUn_JTFH01001739v1_decoy
    chrUn_JTFH01001740v1_decoy
    chrUn_JTFH01001741v1_decoy
    chrUn_JTFH01001742v1_decoy
    chrUn_JTFH01001743v1_decoy
    chrUn_JTFH01001744v1_decoy
    chrUn_JTFH01001745v1_decoy
    chrUn_JTFH01001746v1_decoy
    chrUn_JTFH01001747v1_decoy
    chrUn_JTFH01001748v1_decoy
    chrUn_JTFH01001749v1_decoy
    chrUn_JTFH01001750v1_decoy
    chrUn_JTFH01001751v1_decoy
    chrUn_JTFH01001752v1_decoy
    chrUn_JTFH01001753v1_decoy
    chrUn_JTFH01001754v1_decoy
    chrUn_JTFH01001755v1_decoy
    chrUn_JTFH01001756v1_decoy
    chrUn_JTFH01001757v1_decoy
    chrUn_JTFH01001758v1_decoy
    chrUn_JTFH01001759v1_decoy
    chrUn_JTFH01001760v1_decoy
    chrUn_JTFH01001761v1_decoy
    chrUn_JTFH01001762v1_decoy
    chrUn_JTFH01001763v1_decoy
    chrUn_JTFH01001764v1_decoy
    chrUn_JTFH01001765v1_decoy
    chrUn_JTFH01001766v1_decoy
    chrUn_JTFH01001767v1_decoy
    chrUn_JTFH01001768v1_decoy
    chrUn_JTFH01001769v1_decoy
    chrUn_JTFH01001770v1_decoy
    chrUn_JTFH01001771v1_decoy
    chrUn_JTFH01001772v1_decoy
    chrUn_JTFH01001773v1_decoy
    chrUn_JTFH01001774v1_decoy
    chrUn_JTFH01001775v1_decoy
    chrUn_JTFH01001776v1_decoy
    chrUn_JTFH01001777v1_decoy
    chrUn_JTFH01001778v1_decoy
    chrUn_JTFH01001779v1_decoy
    chrUn_JTFH01001780v1_decoy
    chrUn_JTFH01001781v1_decoy
    chrUn_JTFH01001782v1_decoy
    chrUn_JTFH01001783v1_decoy
    chrUn_JTFH01001784v1_decoy
    chrUn_JTFH01001785v1_decoy
    chrUn_JTFH01001786v1_decoy
    chrUn_JTFH01001787v1_decoy
    chrUn_JTFH01001788v1_decoy
    chrUn_JTFH01001789v1_decoy
    chrUn_JTFH01001790v1_decoy
    chrUn_JTFH01001791v1_decoy
    chrUn_JTFH01001792v1_decoy
    chrUn_JTFH01001793v1_decoy
    chrUn_JTFH01001794v1_decoy
    chrUn_JTFH01001795v1_decoy
    chrUn_JTFH01001796v1_decoy
    chrUn_JTFH01001797v1_decoy
    chrUn_JTFH01001798v1_decoy
    chrUn_JTFH01001799v1_decoy
    chrUn_JTFH01001800v1_decoy
    chrUn_JTFH01001801v1_decoy
    chrUn_JTFH01001802v1_decoy
    chrUn_JTFH01001803v1_decoy
    chrUn_JTFH01001804v1_decoy
    chrUn_JTFH01001805v1_decoy
    chrUn_JTFH01001806v1_decoy
    chrUn_JTFH01001807v1_decoy
    chrUn_JTFH01001808v1_decoy
    chrUn_JTFH01001809v1_decoy
    chrUn_JTFH01001810v1_decoy
    chrUn_JTFH01001811v1_decoy
    chrUn_JTFH01001812v1_decoy
    chrUn_JTFH01001813v1_decoy
    chrUn_JTFH01001814v1_decoy
    chrUn_JTFH01001815v1_decoy
    chrUn_JTFH01001816v1_decoy
    chrUn_JTFH01001817v1_decoy
    chrUn_JTFH01001818v1_decoy
    chrUn_JTFH01001819v1_decoy
    chrUn_JTFH01001820v1_decoy
    chrUn_JTFH01001821v1_decoy
    chrUn_JTFH01001822v1_decoy
    chrUn_JTFH01001823v1_decoy
    chrUn_JTFH01001824v1_decoy
    chrUn_JTFH01001825v1_decoy
    chrUn_JTFH01001826v1_decoy
    chrUn_JTFH01001827v1_decoy
    chrUn_JTFH01001828v1_decoy
    chrUn_JTFH01001829v1_decoy
    chrUn_JTFH01001830v1_decoy
    chrUn_JTFH01001831v1_decoy
    chrUn_JTFH01001832v1_decoy
    chrUn_JTFH01001833v1_decoy
    chrUn_JTFH01001834v1_decoy
    chrUn_JTFH01001835v1_decoy
    chrUn_JTFH01001836v1_decoy
    chrUn_JTFH01001837v1_decoy
    chrUn_JTFH01001838v1_decoy
    chrUn_JTFH01001839v1_decoy
    chrUn_JTFH01001840v1_decoy
    chrUn_JTFH01001841v1_decoy
    chrUn_JTFH01001842v1_decoy
    chrUn_JTFH01001843v1_decoy
    chrUn_JTFH01001844v1_decoy
    chrUn_JTFH01001845v1_decoy
    chrUn_JTFH01001846v1_decoy
    chrUn_JTFH01001847v1_decoy
    chrUn_JTFH01001848v1_decoy
    chrUn_JTFH01001849v1_decoy
    chrUn_JTFH01001850v1_decoy
    chrUn_JTFH01001851v1_decoy
    chrUn_JTFH01001852v1_decoy
    chrUn_JTFH01001853v1_decoy
    chrUn_JTFH01001854v1_decoy
    chrUn_JTFH01001855v1_decoy
    chrUn_JTFH01001856v1_decoy
    chrUn_JTFH01001857v1_decoy
    chrUn_JTFH01001858v1_decoy
    chrUn_JTFH01001859v1_decoy
    chrUn_JTFH01001860v1_decoy
    chrUn_JTFH01001861v1_decoy
    chrUn_JTFH01001862v1_decoy
    chrUn_JTFH01001863v1_decoy
    chrUn_JTFH01001864v1_decoy
    chrUn_JTFH01001865v1_decoy
    chrUn_JTFH01001866v1_decoy
    chrUn_JTFH01001867v1_decoy
    chrUn_JTFH01001868v1_decoy
    chrUn_JTFH01001869v1_decoy
    chrUn_JTFH01001870v1_decoy
    chrUn_JTFH01001871v1_decoy
    chrUn_JTFH01001872v1_decoy
    chrUn_JTFH01001873v1_decoy
    chrUn_JTFH01001874v1_decoy
    chrUn_JTFH01001875v1_decoy
    chrUn_JTFH01001876v1_decoy
    chrUn_JTFH01001877v1_decoy
    chrUn_JTFH01001878v1_decoy
    chrUn_JTFH01001879v1_decoy
    chrUn_JTFH01001880v1_decoy
    chrUn_JTFH01001881v1_decoy
    chrUn_JTFH01001882v1_decoy
    chrUn_JTFH01001883v1_decoy
    chrUn_JTFH01001884v1_decoy
    chrUn_JTFH01001885v1_decoy
    chrUn_JTFH01001886v1_decoy
    chrUn_JTFH01001887v1_decoy
    chrUn_JTFH01001888v1_decoy
    chrUn_JTFH01001889v1_decoy
    chrUn_JTFH01001890v1_decoy
    chrUn_JTFH01001891v1_decoy
    chrUn_JTFH01001892v1_decoy
    chrUn_JTFH01001893v1_decoy
    chrUn_JTFH01001894v1_decoy
    chrUn_JTFH01001895v1_decoy
    chrUn_JTFH01001896v1_decoy
    chrUn_JTFH01001897v1_decoy
    chrUn_JTFH01001898v1_decoy
    chrUn_JTFH01001899v1_decoy
    chrUn_JTFH01001900v1_decoy
    chrUn_JTFH01001901v1_decoy
    chrUn_JTFH01001902v1_decoy
    chrUn_JTFH01001903v1_decoy
    chrUn_JTFH01001904v1_decoy
    chrUn_JTFH01001905v1_decoy
    chrUn_JTFH01001906v1_decoy
    chrUn_JTFH01001907v1_decoy
    chrUn_JTFH01001908v1_decoy
    chrUn_JTFH01001909v1_decoy
    chrUn_JTFH01001910v1_decoy
    chrUn_JTFH01001911v1_decoy
    chrUn_JTFH01001912v1_decoy
    chrUn_JTFH01001913v1_decoy
    chrUn_JTFH01001914v1_decoy
    chrUn_JTFH01001915v1_decoy
    chrUn_JTFH01001916v1_decoy
    chrUn_JTFH01001917v1_decoy
    chrUn_JTFH01001918v1_decoy
    chrUn_JTFH01001919v1_decoy
    chrUn_JTFH01001920v1_decoy
    chrUn_JTFH01001921v1_decoy
    chrUn_JTFH01001922v1_decoy
    chrUn_JTFH01001923v1_decoy
    chrUn_JTFH01001924v1_decoy
    chrUn_JTFH01001925v1_decoy
    chrUn_JTFH01001926v1_decoy
    chrUn_JTFH01001927v1_decoy
    chrUn_JTFH01001928v1_decoy
    chrUn_JTFH01001929v1_decoy
    chrUn_JTFH01001930v1_decoy
    chrUn_JTFH01001931v1_decoy
    chrUn_JTFH01001932v1_decoy
    chrUn_JTFH01001933v1_decoy
    chrUn_JTFH01001934v1_decoy
    chrUn_JTFH01001935v1_decoy
    chrUn_JTFH01001936v1_decoy
    chrUn_JTFH01001937v1_decoy
    chrUn_JTFH01001938v1_decoy
    chrUn_JTFH01001939v1_decoy
    chrUn_JTFH01001940v1_decoy
    chrUn_JTFH01001941v1_decoy
    chrUn_JTFH01001942v1_decoy
    chrUn_JTFH01001943v1_decoy
    chrUn_JTFH01001944v1_decoy
    chrUn_JTFH01001945v1_decoy
    chrUn_JTFH01001946v1_decoy
    chrUn_JTFH01001947v1_decoy
    chrUn_JTFH01001948v1_decoy
    chrUn_JTFH01001949v1_decoy
    chrUn_JTFH01001950v1_decoy
    chrUn_JTFH01001951v1_decoy
    chrUn_JTFH01001952v1_decoy
    chrUn_JTFH01001953v1_decoy
    chrUn_JTFH01001954v1_decoy
    chrUn_JTFH01001955v1_decoy
    chrUn_JTFH01001956v1_decoy
    chrUn_JTFH01001957v1_decoy
    chrUn_JTFH01001958v1_decoy
    chrUn_JTFH01001959v1_decoy
    chrUn_JTFH01001960v1_decoy
    chrUn_JTFH01001961v1_decoy
    chrUn_JTFH01001962v1_decoy
    chrUn_JTFH01001963v1_decoy
    chrUn_JTFH01001964v1_decoy
    chrUn_JTFH01001965v1_decoy
    chrUn_JTFH01001966v1_decoy
    chrUn_JTFH01001967v1_decoy
    chrUn_JTFH01001968v1_decoy
    chrUn_JTFH01001969v1_decoy
    chrUn_JTFH01001970v1_decoy
    chrUn_JTFH01001971v1_decoy
    chrUn_JTFH01001972v1_decoy
    chrUn_JTFH01001973v1_decoy
    chrUn_JTFH01001974v1_decoy
    chrUn_JTFH01001975v1_decoy
    chrUn_JTFH01001976v1_decoy
    chrUn_JTFH01001977v1_decoy
    chrUn_JTFH01001978v1_decoy
    chrUn_JTFH01001979v1_decoy
    chrUn_JTFH01001980v1_decoy
    chrUn_JTFH01001981v1_decoy
    chrUn_JTFH01001982v1_decoy
    chrUn_JTFH01001983v1_decoy
    chrUn_JTFH01001984v1_decoy
    chrUn_JTFH01001985v1_decoy
    chrUn_JTFH01001986v1_decoy
    chrUn_JTFH01001987v1_decoy
    chrUn_JTFH01001988v1_decoy
    chrUn_JTFH01001989v1_decoy
    chrUn_JTFH01001990v1_decoy
    chrUn_JTFH01001991v1_decoy
    chrUn_JTFH01001992v1_decoy
    chrUn_JTFH01001993v1_decoy
    chrUn_JTFH01001994v1_decoy
    chrUn_JTFH01001995v1_decoy
    chrUn_JTFH01001996v1_decoy
    chrUn_JTFH01001997v1_decoy
    chrUn_JTFH01001998v1_decoy
    HLA-A*01:01:01:01
    HLA-A*01:01:01:02N
    HLA-A*01:01:38L
    HLA-A*01:02
    HLA-A*01:03
    HLA-A*01:04N
    HLA-A*01:09
    HLA-A*01:11N
    HLA-A*01:14
    HLA-A*01:16N
    HLA-A*01:20
    HLA-A*02:01:01:01
    HLA-A*02:01:01:02L
    HLA-A*02:01:01:03
    HLA-A*02:01:01:04
    HLA-A*02:02:01
    HLA-A*02:03:01
    HLA-A*02:03:03
    HLA-A*02:05:01
    HLA-A*02:06:01
    HLA-A*02:07:01
    HLA-A*02:10
    HLA-A*02:251
    HLA-A*02:259
    HLA-A*02:264
    HLA-A*02:265
    HLA-A*02:266
    HLA-A*02:269
    HLA-A*02:279
    HLA-A*02:32N
    HLA-A*02:376
    HLA-A*02:43N
    HLA-A*02:455
    HLA-A*02:48
    HLA-A*02:51
    HLA-A*02:533
    HLA-A*02:53N
    HLA-A*02:57
    HLA-A*02:60:01
    HLA-A*02:65
    HLA-A*02:68
    HLA-A*02:77
    HLA-A*02:81
    HLA-A*02:89
    HLA-A*02:95
    HLA-A*03:01:01:01
    HLA-A*03:01:01:02N
    HLA-A*03:01:01:03
    HLA-A*03:02:01
    HLA-A*03:11N
    HLA-A*03:21N
    HLA-A*03:36N
    HLA-A*11:01:01
    HLA-A*11:01:18
    HLA-A*11:02:01
    HLA-A*11:05
    HLA-A*11:110
    HLA-A*11:25
    HLA-A*11:50Q
    HLA-A*11:60
    HLA-A*11:69N
    HLA-A*11:74
    HLA-A*11:75
    HLA-A*11:77
    HLA-A*23:01:01
    HLA-A*23:09
    HLA-A*23:38N
    HLA-A*24:02:01:01
    HLA-A*24:02:01:02L
    HLA-A*24:02:01:03
    HLA-A*24:02:03Q
    HLA-A*24:02:10
    HLA-A*24:03:01
    HLA-A*24:07:01
    HLA-A*24:08
    HLA-A*24:09N
    HLA-A*24:10:01
    HLA-A*24:11N
    HLA-A*24:152
    HLA-A*24:20
    HLA-A*24:215
    HLA-A*24:61
    HLA-A*24:86N
    HLA-A*25:01:01
    HLA-A*26:01:01
    HLA-A*26:11N
    HLA-A*26:15
    HLA-A*26:50
    HLA-A*29:01:01:01
    HLA-A*29:01:01:02N
    HLA-A*29:02:01:01
    HLA-A*29:02:01:02
    HLA-A*29:46
    HLA-A*30:01:01
    HLA-A*30:02:01:01
    HLA-A*30:02:01:02
    HLA-A*30:04:01
    HLA-A*30:89
    HLA-A*31:01:02
    HLA-A*31:01:23
    HLA-A*31:04
    HLA-A*31:14N
    HLA-A*31:46
    HLA-A*32:01:01
    HLA-A*32:06
    HLA-A*33:01:01
    HLA-A*33:03:01
    HLA-A*33:07
    HLA-A*34:01:01
    HLA-A*34:02:01
    HLA-A*36:01
    HLA-A*43:01
    HLA-A*66:01:01
    HLA-A*66:17
    HLA-A*68:01:01:01
    HLA-A*68:01:01:02
    HLA-A*68:01:02:01
    HLA-A*68:01:02:02
    HLA-A*68:02:01:01
    HLA-A*68:02:01:02
    HLA-A*68:02:01:03
    HLA-A*68:02:02
    HLA-A*68:03:01
    HLA-A*68:08:01
    HLA-A*68:113
    HLA-A*68:17
    HLA-A*68:18N
    HLA-A*68:22
    HLA-A*68:71
    HLA-A*69:01
    HLA-A*74:01
    HLA-A*74:02:01:01
    HLA-A*74:02:01:02
    HLA-A*80:01:01:01
    HLA-A*80:01:01:02
    HLA-B*07:02:01
    HLA-B*07:05:01
    HLA-B*07:06
    HLA-B*07:156
    HLA-B*07:33:01
    HLA-B*07:41
    HLA-B*07:44
    HLA-B*07:50
    HLA-B*08:01:01
    HLA-B*08:08N
    HLA-B*08:132
    HLA-B*08:134
    HLA-B*08:19N
    HLA-B*08:20
    HLA-B*08:33
    HLA-B*08:79
    HLA-B*13:01:01
    HLA-B*13:02:01
    HLA-B*13:02:03
    HLA-B*13:02:09
    HLA-B*13:08
    HLA-B*13:15
    HLA-B*13:25
    HLA-B*14:01:01
    HLA-B*14:02:01
    HLA-B*14:07N
    HLA-B*15:01:01:01
    HLA-B*15:01:01:02N
    HLA-B*15:01:01:03
    HLA-B*15:02:01
    HLA-B*15:03:01
    HLA-B*15:04:01
    HLA-B*15:07:01
    HLA-B*15:108
    HLA-B*15:10:01
    HLA-B*15:11:01
    HLA-B*15:13:01
    HLA-B*15:16:01
    HLA-B*15:17:01:01
    HLA-B*15:17:01:02
    HLA-B*15:18:01
    HLA-B*15:220
    HLA-B*15:25:01
    HLA-B*15:27:01
    HLA-B*15:32:01
    HLA-B*15:42
    HLA-B*15:58
    HLA-B*15:66
    HLA-B*15:77
    HLA-B*15:83
    HLA-B*18:01:01:01
    HLA-B*18:01:01:02
    HLA-B*18:02
    HLA-B*18:03
    HLA-B*18:17N
    HLA-B*18:26
    HLA-B*18:94N
    HLA-B*27:04:01
    HLA-B*27:05:02
    HLA-B*27:05:18
    HLA-B*27:06
    HLA-B*27:07:01
    HLA-B*27:131
    HLA-B*27:24
    HLA-B*27:25
    HLA-B*27:32
    HLA-B*35:01:01:01
    HLA-B*35:01:01:02
    HLA-B*35:01:22
    HLA-B*35:02:01
    HLA-B*35:03:01
    HLA-B*35:05:01
    HLA-B*35:08:01
    HLA-B*35:14:02
    HLA-B*35:241
    HLA-B*35:41
    HLA-B*37:01:01
    HLA-B*37:01:05
    HLA-B*38:01:01
    HLA-B*38:02:01
    HLA-B*38:14
    HLA-B*39:01:01:01
    HLA-B*39:01:01:02L
    HLA-B*39:01:01:03
    HLA-B*39:01:03
    HLA-B*39:01:16
    HLA-B*39:01:21
    HLA-B*39:05:01
    HLA-B*39:06:02
    HLA-B*39:10:01
    HLA-B*39:13:02
    HLA-B*39:14
    HLA-B*39:34
    HLA-B*39:38Q
    HLA-B*40:01:01
    HLA-B*40:01:02
    HLA-B*40:02:01
    HLA-B*40:03
    HLA-B*40:06:01:01
    HLA-B*40:06:01:02
    HLA-B*40:10:01
    HLA-B*40:150
    HLA-B*40:40
    HLA-B*40:72:01
    HLA-B*40:79
    HLA-B*41:01:01
    HLA-B*41:02:01
    HLA-B*42:01:01
    HLA-B*42:02
    HLA-B*42:08
    HLA-B*44:02:01:01
    HLA-B*44:02:01:02S
    HLA-B*44:02:01:03
    HLA-B*44:02:17
    HLA-B*44:02:27
    HLA-B*44:03:01
    HLA-B*44:03:02
    HLA-B*44:04
    HLA-B*44:09
    HLA-B*44:138Q
    HLA-B*44:150
    HLA-B*44:23N
    HLA-B*44:26
    HLA-B*44:46
    HLA-B*44:49
    HLA-B*44:56N
    HLA-B*45:01:01
    HLA-B*45:04
    HLA-B*46:01:01
    HLA-B*46:01:05
    HLA-B*47:01:01:01
    HLA-B*47:01:01:02
    HLA-B*48:01:01
    HLA-B*48:03:01
    HLA-B*48:04
    HLA-B*48:08
    HLA-B*49:01:01
    HLA-B*49:32
    HLA-B*50:01:01
    HLA-B*51:01:01
    HLA-B*51:01:02
    HLA-B*51:02:01
    HLA-B*51:07:01
    HLA-B*51:42
    HLA-B*52:01:01:01
    HLA-B*52:01:01:02
    HLA-B*52:01:01:03
    HLA-B*52:01:02
    HLA-B*53:01:01
    HLA-B*53:11
    HLA-B*54:01:01
    HLA-B*54:18
    HLA-B*55:01:01
    HLA-B*55:01:03
    HLA-B*55:02:01
    HLA-B*55:12
    HLA-B*55:24
    HLA-B*55:48
    HLA-B*56:01:01
    HLA-B*56:03
    HLA-B*56:04
    HLA-B*57:01:01
    HLA-B*57:03:01
    HLA-B*57:06
    HLA-B*57:11
    HLA-B*57:29
    HLA-B*58:01:01
    HLA-B*58:31N
    HLA-B*59:01:01:01
    HLA-B*59:01:01:02
    HLA-B*67:01:01
    HLA-B*67:01:02
    HLA-B*67:02
    HLA-B*73:01
    HLA-B*78:01:01
    HLA-B*81:01
    HLA-B*82:02:01
    HLA-C*01:02:01
    HLA-C*01:02:11
    HLA-C*01:02:29
    HLA-C*01:02:30
    HLA-C*01:03
    HLA-C*01:06
    HLA-C*01:08
    HLA-C*01:14
    HLA-C*01:21
    HLA-C*01:30
    HLA-C*01:40
    HLA-C*02:02:02:01
    HLA-C*02:02:02:02
    HLA-C*02:10
    HLA-C*02:11
    HLA-C*02:16:02
    HLA-C*02:69
    HLA-C*02:85
    HLA-C*02:86
    HLA-C*02:87
    HLA-C*03:02:01
    HLA-C*03:02:02:01
    HLA-C*03:02:02:02
    HLA-C*03:02:02:03
    HLA-C*03:03:01
    HLA-C*03:04:01:01
    HLA-C*03:04:01:02
    HLA-C*03:04:02
    HLA-C*03:04:04
    HLA-C*03:05
    HLA-C*03:06
    HLA-C*03:100
    HLA-C*03:13:01
    HLA-C*03:20N
    HLA-C*03:219
    HLA-C*03:261
    HLA-C*03:40:01
    HLA-C*03:41:02
    HLA-C*03:46
    HLA-C*03:61
    HLA-C*04:01:01:01
    HLA-C*04:01:01:02
    HLA-C*04:01:01:03
    HLA-C*04:01:01:04
    HLA-C*04:01:01:05
    HLA-C*04:01:62
    HLA-C*04:03:01
    HLA-C*04:06
    HLA-C*04:09N
    HLA-C*04:128
    HLA-C*04:161
    HLA-C*04:177
    HLA-C*04:70
    HLA-C*04:71
    HLA-C*05:01:01:01
    HLA-C*05:01:01:02
    HLA-C*05:08
    HLA-C*05:09:01
    HLA-C*05:93
    HLA-C*06:02:01:01
    HLA-C*06:02:01:02
    HLA-C*06:02:01:03
    HLA-C*06:23
    HLA-C*06:24
    HLA-C*06:46N
    HLA-C*07:01:01:01
    HLA-C*07:01:01:02
    HLA-C*07:01:02
    HLA-C*07:01:19
    HLA-C*07:01:27
    HLA-C*07:01:45
    HLA-C*07:02:01:01
    HLA-C*07:02:01:02
    HLA-C*07:02:01:03
    HLA-C*07:02:01:04
    HLA-C*07:02:01:05
    HLA-C*07:02:05
    HLA-C*07:02:06
    HLA-C*07:02:64
    HLA-C*07:04:01
    HLA-C*07:04:02
    HLA-C*07:06
    HLA-C*07:149
    HLA-C*07:18
    HLA-C*07:19
    HLA-C*07:26
    HLA-C*07:30
    HLA-C*07:32N
    HLA-C*07:384
    HLA-C*07:385
    HLA-C*07:386
    HLA-C*07:391
    HLA-C*07:392
    HLA-C*07:49
    HLA-C*07:56:02
    HLA-C*07:66
    HLA-C*07:67
    HLA-C*08:01:01
    HLA-C*08:01:03
    HLA-C*08:02:01:01
    HLA-C*08:02:01:02
    HLA-C*08:03:01
    HLA-C*08:04:01
    HLA-C*08:112
    HLA-C*08:20
    HLA-C*08:21
    HLA-C*08:22
    HLA-C*08:24
    HLA-C*08:27
    HLA-C*08:36N
    HLA-C*08:40
    HLA-C*08:41
    HLA-C*08:62
    HLA-C*12:02:02
    HLA-C*12:03:01:01
    HLA-C*12:03:01:02
    HLA-C*12:08
    HLA-C*12:13
    HLA-C*12:19
    HLA-C*12:22
    HLA-C*12:99
    HLA-C*14:02:01
    HLA-C*14:03
    HLA-C*14:21N
    HLA-C*14:23
    HLA-C*15:02:01
    HLA-C*15:05:01
    HLA-C*15:05:02
    HLA-C*15:13
    HLA-C*15:16
    HLA-C*15:17
    HLA-C*15:96Q
    HLA-C*16:01:01
    HLA-C*16:02:01
    HLA-C*16:04:01
    HLA-C*17:01:01:01
    HLA-C*17:01:01:02
    HLA-C*17:01:01:03
    HLA-C*17:03
    HLA-C*18:01
    HLA-DQA1*01:01:02
    HLA-DQA1*01:02:01:01
    HLA-DQA1*01:02:01:02
    HLA-DQA1*01:02:01:03
    HLA-DQA1*01:02:01:04
    HLA-DQA1*01:03:01:01
    HLA-DQA1*01:03:01:02
    HLA-DQA1*01:04:01:01
    HLA-DQA1*01:04:01:02
    HLA-DQA1*01:05:01
    HLA-DQA1*01:07
    HLA-DQA1*01:10
    HLA-DQA1*01:11
    HLA-DQA1*02:01
    HLA-DQA1*03:01:01
    HLA-DQA1*03:02
    HLA-DQA1*03:03:01
    HLA-DQA1*04:01:02:01
    HLA-DQA1*04:01:02:02
    HLA-DQA1*04:02
    HLA-DQA1*05:01:01:01
    HLA-DQA1*05:01:01:02
    HLA-DQA1*05:03
    HLA-DQA1*05:05:01:01
    HLA-DQA1*05:05:01:02
    HLA-DQA1*05:05:01:03
    HLA-DQA1*05:11
    HLA-DQA1*06:01:01
    HLA-DQB1*02:01:01
    HLA-DQB1*02:02:01
    HLA-DQB1*03:01:01:01
    HLA-DQB1*03:01:01:02
    HLA-DQB1*03:01:01:03
    HLA-DQB1*03:02:01
    HLA-DQB1*03:03:02:01
    HLA-DQB1*03:03:02:02
    HLA-DQB1*03:03:02:03
    HLA-DQB1*03:05:01
    HLA-DQB1*05:01:01:01
    HLA-DQB1*05:01:01:02
    HLA-DQB1*05:03:01:01
    HLA-DQB1*05:03:01:02
    HLA-DQB1*06:01:01
    HLA-DQB1*06:02:01
    HLA-DQB1*06:03:01
    HLA-DQB1*06:09:01
    HLA-DRB1*01:01:01
    HLA-DRB1*01:02:01
    HLA-DRB1*03:01:01:01
    HLA-DRB1*03:01:01:02
    HLA-DRB1*04:03:01
    HLA-DRB1*07:01:01:01
    HLA-DRB1*07:01:01:02
    HLA-DRB1*08:03:02
    HLA-DRB1*09:21
    HLA-DRB1*10:01:01
    HLA-DRB1*11:01:01
    HLA-DRB1*11:01:02
    HLA-DRB1*11:04:01
    HLA-DRB1*12:01:01
    HLA-DRB1*12:17
    HLA-DRB1*13:01:01
    HLA-DRB1*13:02:01
    HLA-DRB1*14:05:01
    HLA-DRB1*14:54:01
    HLA-DRB1*15:01:01:01
    HLA-DRB1*15:01:01:02
    HLA-DRB1*15:01:01:03
    HLA-DRB1*15:01:01:04
    HLA-DRB1*15:02:01
    HLA-DRB1*15:03:01:01
    HLA-DRB1*15:03:01:02
    HLA-DRB1*16:02:01
--}
