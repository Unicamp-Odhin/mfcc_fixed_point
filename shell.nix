{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  packages = with pkgs; [
    zlib
    gcc
    gnumake
  ];

  C_INCLUDE_PATH = "${pkgs.zlib}/include";
  LIBRARY_PATH = "${pkgs.zlib}/lib";
}
