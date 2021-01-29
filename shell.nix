# shell.nix
let
  sources = import ./nix/sources.nix;
  rustPlatform = import ./nix/rust.nix { inherit sources; };
  pkgs = import sources.nixpkgs {};
in
pkgs.mkShell {
  buildInputs = [
    rustPlatform.rust.cargo
    pkgs.cargo-edit
  ];
}
