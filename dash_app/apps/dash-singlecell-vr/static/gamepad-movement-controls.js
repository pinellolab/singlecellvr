AFRAME.registerComponent('gamepad-movement-controls', {
    init: function () {
        let controlsEl = document.querySelector('[button-controls]');
        this.moving = false;
        this.onbuttondown = function(e) {
          this.toggleMovement();
        }
        this.onbuttonup = function(e) {
          this.toggleMovement();
        }
        controlsEl.addEventListener('buttondown', this.onbuttondown.bind(this));
        controlsEl.addEventListener('buttonup', this.onbuttonup.bind(this));
    },
    tick: function (time, timeDelta) {
        if (this.moving) {
            // Duplicated in index.js, figure out how to unify these
            let direction = new THREE.Vector3();
            const camera = AFRAME.scenes[0].camera;
            camera.getWorldDirection( direction );
            direction.multiplyScalar(.05);
            const cameraEl = document.getElementById('rig');
            var pos = cameraEl.getAttribute("position");
            pos.x += direction.x
            pos.y += direction.y
            pos.z += direction.z
            const mapPlayer = document.getElementById('mapPlayer').object3D;
            mapPlayer.position.set((pos.x + direction.x)  * .01, (pos.y + direction.y) * .01, (pos.z + direction.z) * .01);
        }
    },
    toggleMovement: function () {
      this.moving = !this.moving;
    }
});