import os
import psutil

from ctypes import cdll


from OCC.Core.gp import gp_Pnt, gp_Vec

pnt = gp_Pnt(0, 0, 0)

# 現在のプロセスIDを取得
pid = os.getpid()

# psutilでプロセス情報を取得
process = psutil.Process(pid)

# ロードされているモジュール（DLL）をリストアップ
print("Loaded DLLs:")
for dll in process.memory_maps():
    if dll.path and "TK" in dll.path:  # OpenCASCADE関連のDLLをフィルタリング
        # DLLをロード
        dll_c = cdll.LoadLibrary(dll.path)
        print(f"Loaded DLL: {dll.path}")
        print(f"          : {dll_c}")

# VS Developer Command Promptを使用して、以下のコマンドを実行してDLLの依存関係を確認
# dumpbin /DEPENDENTS <DLLファイル名>
